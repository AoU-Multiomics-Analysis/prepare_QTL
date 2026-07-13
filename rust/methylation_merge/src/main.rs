use clap::Parser;
use csv::{Reader, ReaderBuilder, StringRecord, Writer, WriterBuilder};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read};
use std::path::{Path, PathBuf};
use thiserror::Error;

const MAX_OPEN_STREAMS: usize = 128;
const GZIP_BUFFER: usize = 1024 * 1024;

#[derive(Parser, Debug)]
#[command(about = "Bounded-memory cohort merge for one methylation chromosome")]
struct Args {
    #[arg(long)]
    all_call_list: PathBuf,
    #[arg(long)]
    sample_qc: PathBuf,
    #[arg(long)]
    cohort_samples: PathBuf,
    #[arg(long)]
    total_samples: usize,
    #[arg(long)]
    chromosome: String,
    #[arg(long)]
    output_prefix: String,
    #[arg(long, default_value_t = 0.95)]
    min_sample_fraction: f64,
    #[arg(long, default_value_t = 0)]
    min_samples: usize,
    #[arg(long, default_value_t = 0.003)]
    min_methylation_mad: f64,
    #[arg(long, default_value = "mod_score")]
    value_column: String,
    #[arg(long, default_value_t = 0.01)]
    value_multiplier: f64,
    #[arg(long, default_value_t = 1000)]
    progress_every_sites: usize,
}

#[derive(Error, Debug)]
enum MergeError {
    #[error("{0}")]
    Message(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Csv(#[from] csv::Error),
}
type Result<T> = std::result::Result<T, MergeError>;

#[derive(Clone)]
struct Schema {
    header: StringRecord,
    chrom: usize,
    begin: usize,
    end: usize,
    site_key: usize,
    sample_id: usize,
    cov: usize,
    meets_min_coverage: usize,
    per_sample_qc_pass: usize,
    value: usize,
}

#[derive(Clone, Eq, PartialEq)]
struct SortKey {
    begin: i64,
    end: i64,
}
impl Ord for SortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.begin.cmp(&other.begin).then(self.end.cmp(&other.end))
    }
}
impl PartialOrd for SortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct CallStream {
    reader: Reader<Box<dyn Read>>,
    schema: Schema,
    previous: Option<SortKey>,
    path: PathBuf,
}
#[derive(Eq, PartialEq)]
struct HeapEntry {
    key: SortKey,
    stream: usize,
    record: StringRecord,
}
impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .key
            .cmp(&self.key)
            .then(other.stream.cmp(&self.stream))
    }
}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
struct KWayMerge {
    streams: Vec<CallStream>,
    heap: BinaryHeap<HeapEntry>,
}

#[derive(Clone)]
struct Call {
    record: StringRecord,
    sample_index: usize,
    cov: Option<f64>,
    value: Option<f64>,
    min_cov: bool,
    pass: bool,
    normalized_log_cov: Option<f64>,
}

fn open_read(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path)?;
    if path
        .extension()
        .is_some_and(|x| x.eq_ignore_ascii_case("gz"))
    {
        Ok(Box::new(MultiGzDecoder::new(BufReader::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}
fn open_writer(path: &Path) -> Result<Writer<BufWriter<GzEncoder<File>>>> {
    let encoder = GzEncoder::new(File::create(path)?, Compression::default());
    Ok(WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(BufWriter::with_capacity(GZIP_BUFFER, encoder)))
}
fn column(header: &StringRecord, name: &str) -> Result<usize> {
    header.iter().position(|x| x == name).ok_or_else(|| {
        MergeError::Message(format!(
            "Input call tables are missing required column '{name}'"
        ))
    })
}
fn schema_from_header(header: StringRecord, value_column: &str) -> Result<Schema> {
    Ok(Schema {
        chrom: column(&header, "#chrom")?,
        begin: column(&header, "begin")?,
        end: column(&header, "end")?,
        site_key: column(&header, "site_key")?,
        sample_id: column(&header, "sample_id")?,
        cov: column(&header, "cov")?,
        meets_min_coverage: column(&header, "meets_min_coverage")?,
        per_sample_qc_pass: column(&header, "per_sample_qc_pass")?,
        value: column(&header, value_column)?,
        header,
    })
}
fn key(record: &StringRecord, schema: &Schema, path: &Path) -> Result<SortKey> {
    let parse = |index: usize, label: &str| {
        record
            .get(index)
            .unwrap_or_default()
            .parse::<i64>()
            .map_err(|_| MergeError::Message(format!("Invalid {label} in {}", path.display())))
    };
    Ok(SortKey {
        begin: parse(schema.begin, "begin")?,
        end: parse(schema.end, "end")?,
    })
}
impl CallStream {
    fn open(path: PathBuf, schema: &Schema) -> Result<Self> {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .from_reader(open_read(&path)?);
        let header = reader.headers()?.clone();
        if header != schema.header {
            return Err(MergeError::Message(format!(
                "All-call files do not use the same schema: {}",
                path.display()
            )));
        }
        Ok(Self {
            reader,
            schema: schema.clone(),
            previous: None,
            path,
        })
    }
    fn next_record(&mut self) -> Result<Option<(SortKey, StringRecord)>> {
        let mut record = StringRecord::new();
        if !self.reader.read_record(&mut record)? {
            return Ok(None);
        }
        let current = key(&record, &self.schema, &self.path)?;
        if let Some(previous) = &self.previous {
            if current < *previous {
                return Err(MergeError::Message(format!("{} is not sorted by begin and end. The Rust cohort merger requires one coordinate-sorted file per sample.", self.path.display())));
            }
        }
        self.previous = Some(current.clone());
        Ok(Some((current, record)))
    }
}
impl KWayMerge {
    fn new(paths: &[PathBuf], schema: &Schema) -> Result<Self> {
        if paths.len() > MAX_OPEN_STREAMS {
            return Err(MergeError::Message(format!(
                "Internal merge received {} streams; expected at most {MAX_OPEN_STREAMS}",
                paths.len()
            )));
        }
        let mut merger = Self {
            streams: paths
                .iter()
                .cloned()
                .map(|path| CallStream::open(path, schema))
                .collect::<Result<Vec<_>>>()?,
            heap: BinaryHeap::new(),
        };
        for index in 0..merger.streams.len() {
            merger.push_next(index)?;
        }
        Ok(merger)
    }
    fn push_next(&mut self, index: usize) -> Result<()> {
        if let Some((key, record)) = self.streams[index].next_record()? {
            self.heap.push(HeapEntry {
                key,
                stream: index,
                record,
            });
        }
        Ok(())
    }
    fn next(&mut self) -> Result<Option<(SortKey, StringRecord)>> {
        let Some(entry) = self.heap.pop() else {
            return Ok(None);
        };
        self.push_next(entry.stream)?;
        Ok(Some((entry.key, entry.record)))
    }
}

fn path_list(path: &Path) -> Result<Vec<PathBuf>> {
    let content = fs::read_to_string(path)?;
    let root = path.parent().unwrap_or_else(|| Path::new("."));
    let paths: Vec<_> = content
        .lines()
        .filter(|x| !x.trim().is_empty())
        .map(|x| {
            let p = PathBuf::from(x.trim());
            if p.is_absolute() {
                p
            } else {
                root.join(p)
            }
        })
        .collect();
    if paths.is_empty() {
        return Err(MergeError::Message(format!(
            "All-call file list is empty: {}",
            path.display()
        )));
    }
    for path in &paths {
        if !path.exists() {
            return Err(MergeError::Message(format!(
                "All-call file does not exist: {}",
                path.display()
            )));
        }
    }
    Ok(paths)
}
fn first_schema(paths: &[PathBuf], value_column: &str) -> Result<Schema> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(open_read(&paths[0])?);
    schema_from_header(reader.headers()?.clone(), value_column)
}
fn merge_run(paths: &[PathBuf], schema: &Schema, output: &Path) -> Result<()> {
    let mut writer = open_writer(output)?;
    writer.write_record(&schema.header)?;
    let mut merge = KWayMerge::new(paths, schema)?;
    while let Some((_key, record)) = merge.next()? {
        writer.write_record(&record)?;
    }
    writer.flush()?;
    Ok(())
}
fn reduce_paths(paths: Vec<PathBuf>, schema: &Schema, temp_dir: &Path) -> Result<Vec<PathBuf>> {
    if paths.len() <= MAX_OPEN_STREAMS {
        return Ok(paths);
    }
    fs::create_dir_all(temp_dir)?;
    let mut runs = Vec::new();
    for (index, chunk) in paths.chunks(MAX_OPEN_STREAMS).enumerate() {
        let run = temp_dir.join(format!("sorted_run_{index:05}.tsv.gz"));
        eprintln!(
            "Creating sorted intermediate run {}/{} from {} input files",
            index + 1,
            paths.len().div_ceil(MAX_OPEN_STREAMS),
            chunk.len()
        );
        merge_run(chunk, schema, &run)?;
        runs.push(run);
    }
    Ok(runs)
}

fn parse_number(value: &str, label: &str) -> Result<Option<f64>> {
    if value.is_empty() || value == "NA" || value == "." {
        return Ok(None);
    }
    value
        .parse::<f64>()
        .map(Some)
        .map_err(|_| MergeError::Message(format!("Column '{label}' must be numeric")))
}
fn parse_bool(value: &str, label: &str) -> Result<bool> {
    match value {
        "TRUE" | "true" | "T" => Ok(true),
        "FALSE" | "false" | "F" => Ok(false),
        _ => Err(MergeError::Message(format!(
            "All-call files contain non-logical {label} values"
        ))),
    }
}
fn number(value: Option<f64>) -> String {
    match value {
        Some(x) if x.is_finite() => x.to_string(),
        _ => "NA".to_owned(),
    }
}
fn bool_text(value: bool) -> &'static str {
    if value {
        "TRUE"
    } else {
        "FALSE"
    }
}
fn finite(values: impl Iterator<Item = Option<f64>>) -> Vec<f64> {
    values.filter_map(|x| x.filter(|v| v.is_finite())).collect()
}
fn mean(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        None
    } else {
        Some(values.iter().sum::<f64>() / values.len() as f64)
    }
}
struct NumericSummary {
    mean: Option<f64>,
    sd: Option<f64>,
    cv: Option<f64>,
}
fn summarize(values: &[f64]) -> NumericSummary {
    let average = mean(values);
    let standard_deviation = if values.len() < 2 {
        None
    } else {
        average.map(|avg| {
            (values.iter().map(|x| (x - avg).powi(2)).sum::<f64>() / (values.len() - 1) as f64)
                .sqrt()
        })
    };
    let coefficient_of_variation = match (average, standard_deviation) {
        (Some(avg), Some(std)) if avg != 0.0 => Some(std / avg),
        _ => None,
    };
    NumericSummary {
        mean: average,
        sd: standard_deviation,
        cv: coefficient_of_variation,
    }
}
fn quantile(values: &mut [f64], probability: f64) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let index = (values.len() - 1) as f64 * probability;
    let low = index.floor() as usize;
    let high = index.ceil() as usize;
    Some(values[low] + (index - low as f64) * (values[high] - values[low]))
}
fn median(values: &[f64]) -> Option<f64> {
    quantile(&mut values.to_vec(), 0.5)
}
fn mad(values: &[f64]) -> Option<f64> {
    let med = median(values)?;
    let deviations: Vec<f64> = values.iter().map(|x| (x - med).abs()).collect();
    median(&deviations)
}
fn ranks(values: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|a, b| values[*a].total_cmp(&values[*b]));
    let mut result = vec![0.0; values.len()];
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && values[order[end]] == values[order[start]] {
            end += 1;
        }
        let rank = (start + 1 + end) as f64 / 2.0;
        for index in &order[start..end] {
            result[*index] = rank;
        }
        start = end;
    }
    result
}
fn spearman(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.len() < 3 || left.len() != right.len() {
        return None;
    }
    if left.iter().all(|x| *x == left[0]) || right.iter().all(|x| *x == right[0]) {
        return None;
    }
    let x = ranks(left);
    let y = ranks(right);
    let xm = mean(&x)?;
    let ym = mean(&y)?;
    let numerator: f64 = x.iter().zip(&y).map(|(a, b)| (a - xm) * (b - ym)).sum();
    let xd: f64 = x.iter().map(|a| (a - xm).powi(2)).sum();
    let yd: f64 = y.iter().map(|b| (b - ym).powi(2)).sum();
    if xd == 0.0 || yd == 0.0 {
        None
    } else {
        Some(numerator / (xd * yd).sqrt())
    }
}
fn inverse_normal(values: &[f64]) -> Vec<f64> {
    let ranks = ranks(values);
    ranks
        .iter()
        .map(|rank| normal_quantile((rank - 0.5) / values.len() as f64))
        .collect()
}
fn normal_quantile(p: f64) -> f64 {
    let a: [f64; 6] = [
        -39.69683028665376,
        220.9460984245205,
        -275.9285104469687,
        138.357751867269,
        -30.66479806614716,
        2.506628277459239,
    ];
    let b: [f64; 5] = [
        -54.47609879822406,
        161.5858368580409,
        -155.6989798598866,
        66.80131188771972,
        -13.28068155288572,
    ];
    let c: [f64; 6] = [
        -0.007784894002430293,
        -0.3223964580411365,
        -2.400758277161838,
        -2.549732539343734,
        4.374664141464968,
        2.938163982698783,
    ];
    let d: [f64; 4] = [
        0.007784695709041462,
        0.3224671290700398,
        2.445134137142996,
        3.754408661907416,
    ];
    if p < 0.02425 {
        let q = (-2.0 * p.ln()).sqrt();
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if p > 0.97575 {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    let q = p - 0.5;
    let r = q * q;
    (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
        / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
}

fn read_cohort(path: &Path, total: usize) -> Result<(Vec<String>, HashMap<String, usize>)> {
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let header = reader.headers()?.clone();
    if header.len() != 1 || header.get(0) != Some("sample_id") {
        return Err(MergeError::Message(
            "--CohortSamples must contain exactly one column named sample_id".to_owned(),
        ));
    }
    let mut samples = Vec::new();
    let mut map = HashMap::new();
    for row in reader.records() {
        let id = row?.get(0).unwrap_or_default().to_owned();
        if id.is_empty() || map.insert(id.clone(), samples.len()).is_some() {
            return Err(MergeError::Message(
                "--CohortSamples contains an empty or duplicate sample_id".to_owned(),
            ));
        }
        samples.push(id);
    }
    if samples.len() != total {
        return Err(MergeError::Message(format!(
            "--CohortSamples contains {} samples, but --TotalSamples is {total}",
            samples.len()
        )));
    }
    Ok((samples, map))
}
fn read_baselines(path: &Path, cohort: &HashMap<String, usize>, total: usize) -> Result<Vec<f64>> {
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let header = reader.headers()?.clone();
    let sid = column(&header, "sample_id")?;
    let median_cov = column(&header, "median_cov")?;
    let mut output = vec![f64::NAN; total];
    for row in reader.records() {
        let row = row?;
        let id = row.get(sid).unwrap_or_default();
        let Some(&index) = cohort.get(id) else {
            return Err(MergeError::Message(format!(
                "Sample-QC contains a sample absent from --CohortSamples: {id}"
            )));
        };
        if output[index].is_finite() {
            return Err(MergeError::Message(
                "Sample-QC contains duplicate sample_id values".to_owned(),
            ));
        }
        let value = parse_number(row.get(median_cov).unwrap_or_default(), "median_cov")?
            .filter(|x| *x > 0.0)
            .ok_or_else(|| {
                MergeError::Message(
                    "Sample-QC contains a missing or non-positive median_cov value".to_owned(),
                )
            })?;
        output[index] = value.ln_1p();
    }
    if output.iter().any(|x| !x.is_finite()) {
        return Err(MergeError::Message(
            "Sample-QC does not contain exactly the cohort sample IDs".to_owned(),
        ));
    }
    Ok(output)
}

struct Outputs {
    long: Writer<BufWriter<GzEncoder<File>>>,
    qc: Writer<BufWriter<GzEncoder<File>>>,
    metadata: Writer<BufWriter<GzEncoder<File>>>,
    raw: Writer<BufWriter<GzEncoder<File>>>,
    int: Writer<BufWriter<GzEncoder<File>>>,
}
struct SiteScratch {
    seen_epoch: Vec<u32>,
    epoch: u32,
}
impl SiteScratch {
    fn new(total_samples: usize) -> Self {
        Self {
            seen_epoch: vec![0; total_samples],
            epoch: 0,
        }
    }
    fn begin_site(&mut self) -> u32 {
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.seen_epoch.fill(0);
            self.epoch = 1;
        }
        self.epoch
    }
    fn mark_seen(&mut self, sample_index: usize, epoch: u32) -> bool {
        if self.seen_epoch[sample_index] == epoch {
            false
        } else {
            self.seen_epoch[sample_index] = epoch;
            true
        }
    }
}
fn write_bed_row(
    writer: &mut Writer<BufWriter<GzEncoder<File>>>,
    chrom: &str,
    key: &SortKey,
    site_key: &str,
    values: &[f64],
) -> Result<()> {
    let begin = key.begin.to_string();
    let end = key.end.to_string();
    writer.write_field(chrom)?;
    writer.write_field(&begin)?;
    writer.write_field(&end)?;
    writer.write_field(site_key)?;
    for value in values {
        writer.write_field(number(Some(*value)))?;
    }
    writer.write_record(None::<&[u8]>)?;
    Ok(())
}
fn outputs(prefix: &str, schema: &Schema, samples: &[String]) -> Result<Outputs> {
    let parent = Path::new(prefix).parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent)?;
    let mut long = open_writer(&PathBuf::from(format!(
        "{prefix}.methylation.filtered.long.tsv.gz"
    )))?;
    let mut long_header = vec!["sample_id".to_owned()];
    long_header.extend(
        schema
            .header
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != schema.sample_id)
            .map(|(_, value)| value.to_owned()),
    );
    long_header.push("methylation_value_for_qtl".to_owned());
    long.write_record(long_header)?;
    let mut qc = open_writer(&PathBuf::from(format!(
        "{prefix}.methylation.site_qc.tsv.gz"
    )))?;
    qc.write_record([
        "#chrom",
        "begin",
        "end",
        "site_key",
        "n_samples_passing",
        "fraction_samples_passing",
        "median_cov_passing",
        "min_cov_passing",
        "max_cov_passing",
        "n_samples_required",
        "keep_site",
    ])?;
    let mut metadata = open_writer(&PathBuf::from(format!(
        "{prefix}.methylation.site_metadata.tsv.gz"
    )))?;
    metadata.write_record([
        "#chrom",
        "begin",
        "end",
        "site_key",
        "n_samples_observed",
        "fraction_samples_observed",
        "mean_cov_all_calls",
        "sd_cov_all_calls",
        "cv_cov_all_calls",
        "mean_methylation_all_calls",
        "sd_methylation_all_calls",
        "cv_methylation_all_calls",
        "n_samples_min_coverage",
        "fraction_samples_min_coverage",
        "n_samples_passing_per_sample_qc",
        "fraction_samples_passing_per_sample_qc",
        "mean_cov_passing_per_sample_qc",
        "sd_cov_passing_per_sample_qc",
        "cv_cov_passing_per_sample_qc",
        "median_cov_passing_per_sample_qc",
        "min_cov_passing_per_sample_qc",
        "max_cov_passing_per_sample_qc",
        "mean_methylation_passing_per_sample_qc",
        "sd_methylation_passing_per_sample_qc",
        "cv_methylation_passing_per_sample_qc",
        "methylation_mad_passing_per_sample_qc",
        "n_samples_coverage_methylation_correlation",
        "coverage_methylation_spearman_rho",
        "n_samples_required",
        "pass_minimum_coverage_filter",
        "pass_sample_presence_filter",
        "pass_methylation_mad_filter",
        "has_missing_or_low_coverage",
        "has_extreme_coverage_loss",
        "keep_site",
        "failure_reason",
        "n_samples_imputed_in_qtl_bed",
    ])?;
    let mut raw = open_writer(&PathBuf::from(format!("{prefix}.methylation.raw.bed.gz")))?;
    let mut int = open_writer(&PathBuf::from(format!("{prefix}.methylation.INT.bed.gz")))?;
    let mut bed_header = vec![
        "#chr".to_owned(),
        "start".to_owned(),
        "end".to_owned(),
        "phenotype_id".to_owned(),
    ];
    bed_header.extend(samples.iter().cloned());
    raw.write_record(&bed_header)?;
    int.write_record(&bed_header)?;
    Ok(Outputs {
        long,
        qc,
        metadata,
        raw,
        int,
    })
}

fn process_site(
    calls: &[Call],
    key: &SortKey,
    schema: &Schema,
    total: usize,
    required: usize,
    min_mad: f64,
    outputs: &mut Outputs,
    scratch: &mut SiteScratch,
) -> Result<bool> {
    let first = calls.first().ok_or_else(|| {
        MergeError::Message("Cannot summarize an empty methylation site".to_owned())
    })?;
    let site_chrom = first.record.get(schema.chrom).unwrap_or_default();
    let site_key = first.record.get(schema.site_key).unwrap_or_default();
    let epoch = scratch.begin_site();
    let mut n_observed = 0usize;
    let mut raw_values = vec![None; total];
    let mut cov_all = Vec::new();
    let mut methyl_all = Vec::new();
    let mut cov_pass = Vec::new();
    let mut methyl_pass = Vec::new();
    let mut corr_methyl = Vec::new();
    let mut corr_cov = Vec::new();
    let mut n_min = 0usize;
    let mut n_pass = 0usize;
    for call in calls {
        let id = call.record.get(schema.sample_id).unwrap_or_default();
        if !scratch.mark_seen(call.sample_index, epoch) {
            return Err(MergeError::Message(format!(
                "Found duplicated sample/site calls for {id} at {}",
                site_key
            )));
        }
        n_observed += 1;
        if let Some(x) = call.cov.filter(|x| x.is_finite()) {
            cov_all.push(x);
        }
        if let Some(x) = call.value.filter(|x| x.is_finite()) {
            methyl_all.push(x);
        }
        if call.min_cov {
            n_min += 1;
        }
        if call.pass {
            n_pass += 1;
            if let Some(x) = call.cov.filter(|x| x.is_finite()) {
                cov_pass.push(x);
            }
            if let Some(x) = call.value.filter(|x| x.is_finite()) {
                methyl_pass.push(x);
            }
            raw_values[call.sample_index] = call.value.filter(|x| x.is_finite());
            if let (Some(m), Some(c)) = (
                call.value.filter(|x| x.is_finite()),
                call.normalized_log_cov,
            ) {
                corr_methyl.push(m);
                corr_cov.push(c);
            }
        }
    }
    let pass_min = n_min >= required;
    let pass_presence = n_pass >= required;
    let cov_all_summary = summarize(&cov_all);
    let methyl_all_summary = summarize(&methyl_all);
    let cov_pass_summary = summarize(&cov_pass);
    let methyl_pass_summary = summarize(&methyl_pass);
    let median_cov_pass = median(&cov_pass);
    let methyl_mad = mad(&methyl_pass);
    let pass_mad = methyl_mad.is_some_and(|x| x >= min_mad);
    let keep = pass_presence && pass_mad;
    let failure = if !pass_min {
        "Insufficient minimum coverage"
    } else if !pass_presence {
        "Extreme coverage exclusion"
    } else if !pass_mad {
        "Low methylation MAD"
    } else {
        "Pass all cohort filters"
    };
    let mut imputed = 0usize;
    if keep {
        let observed = finite(raw_values.iter().copied());
        let average = mean(&observed).ok_or_else(|| {
            MergeError::Message(format!(
                "Cannot impute retained QTL feature with no observed methylation values: {}",
                site_key
            ))
        })?;
        for value in &mut raw_values {
            if value.is_none() {
                *value = Some(average);
                imputed += 1;
            }
        }
        for call in calls.iter().filter(|call| call.pass) {
            outputs
                .long
                .write_field(call.record.get(schema.sample_id).unwrap_or_default())?;
            for (index, value) in call.record.iter().enumerate() {
                if index != schema.sample_id {
                    outputs.long.write_field(value)?;
                }
            }
            outputs.long.write_field(number(call.value))?;
            outputs.long.write_record(None::<&[u8]>)?;
        }
        let values: Vec<f64> = raw_values.iter().map(|x| x.unwrap()).collect();
        let ints = inverse_normal(&values);
        write_bed_row(&mut outputs.raw, site_chrom, key, site_key, &values)?;
        write_bed_row(&mut outputs.int, site_chrom, key, site_key, &ints)?;
    }
    outputs.qc.write_record([
        site_chrom.to_owned(),
        key.begin.to_string(),
        key.end.to_string(),
        site_key.to_owned(),
        n_pass.to_string(),
        (n_pass as f64 / total as f64).to_string(),
        number(median_cov_pass),
        number(cov_pass.iter().copied().reduce(f64::min)),
        number(cov_pass.iter().copied().reduce(f64::max)),
        required.to_string(),
        bool_text(keep).to_owned(),
    ])?;
    outputs.metadata.write_record([
        site_chrom.to_owned(),
        key.begin.to_string(),
        key.end.to_string(),
        site_key.to_owned(),
        n_observed.to_string(),
        (n_observed as f64 / total as f64).to_string(),
        number(cov_all_summary.mean),
        number(cov_all_summary.sd),
        number(cov_all_summary.cv),
        number(methyl_all_summary.mean),
        number(methyl_all_summary.sd),
        number(methyl_all_summary.cv),
        n_min.to_string(),
        (n_min as f64 / total as f64).to_string(),
        n_pass.to_string(),
        (n_pass as f64 / total as f64).to_string(),
        number(cov_pass_summary.mean),
        number(cov_pass_summary.sd),
        number(cov_pass_summary.cv),
        number(median_cov_pass),
        number(cov_pass.iter().copied().reduce(f64::min)),
        number(cov_pass.iter().copied().reduce(f64::max)),
        number(methyl_pass_summary.mean),
        number(methyl_pass_summary.sd),
        number(methyl_pass_summary.cv),
        number(methyl_mad),
        corr_methyl.len().to_string(),
        number(spearman(&corr_methyl, &corr_cov)),
        required.to_string(),
        bool_text(pass_min).to_owned(),
        bool_text(pass_presence).to_owned(),
        bool_text(pass_mad).to_owned(),
        bool_text(n_min < total).to_owned(),
        bool_text(n_pass < n_min).to_owned(),
        bool_text(keep).to_owned(),
        failure.to_owned(),
        imputed.to_string(),
    ])?;
    Ok(keep)
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.total_samples == 0
        || !(0.0 < args.min_sample_fraction && args.min_sample_fraction <= 1.0)
        || args.min_methylation_mad < 0.0
        || args.value_multiplier <= 0.0
        || args.progress_every_sites == 0
    {
        return Err(MergeError::Message(
            "Invalid cohort merge thresholds".to_owned(),
        ));
    }
    let input_paths = path_list(&args.all_call_list)?;
    let schema = first_schema(&input_paths, &args.value_column)?;
    let (samples, sample_index) = read_cohort(&args.cohort_samples, args.total_samples)?;
    let baselines = read_baselines(&args.sample_qc, &sample_index, args.total_samples)?;
    let required = ((args.total_samples as f64 * args.min_sample_fraction).ceil() as usize)
        .max(args.min_samples);
    eprintln!(
        "A site must pass per-sample QC in at least {required} of {} samples",
        args.total_samples
    );
    let temp_dir = PathBuf::from(format!("{}.methylation_merge_tmp", args.output_prefix));
    let paths = reduce_paths(input_paths, &schema, &temp_dir)?;
    let mut merge = KWayMerge::new(&paths, &schema)?;
    let mut outputs = outputs(&args.output_prefix, &schema, &samples)?;
    let mut current: Vec<Call> = Vec::new();
    let mut scratch = SiteScratch::new(args.total_samples);
    let mut current_site: Option<SortKey> = None;
    let mut processed = 0usize;
    let mut kept = 0usize;
    while let Some((sort_key, record)) = merge.next()? {
        if record.get(schema.chrom).unwrap_or_default() != args.chromosome {
            return Err(MergeError::Message(format!(
                "All-call files contain a chromosome other than {}",
                args.chromosome
            )));
        }
        let site = sort_key;
        if current_site.as_ref().is_some_and(|x| x != &site) {
            if process_site(
                &current,
                current_site.as_ref().unwrap(),
                &schema,
                args.total_samples,
                required,
                args.min_methylation_mad,
                &mut outputs,
                &mut scratch,
            )? {
                kept += 1;
            }
            processed += 1;
            if processed % args.progress_every_sites == 0 {
                eprintln!(
                    "Cohort metric progress for {}: {processed} sites processed",
                    args.chromosome
                );
            }
            current.clear();
        }
        let sample_id = record.get(schema.sample_id).unwrap_or_default();
        let index = *sample_index.get(sample_id).ok_or_else(|| {
            MergeError::Message(format!(
                "All-call files contain sample absent from --CohortSamples: {sample_id}"
            ))
        })?;
        let cov = parse_number(record.get(schema.cov).unwrap_or_default(), "cov")?;
        let value = parse_number(
            record.get(schema.value).unwrap_or_default(),
            &args.value_column,
        )?
        .map(|x| x * args.value_multiplier);
        let min_cov = parse_bool(
            record.get(schema.meets_min_coverage).unwrap_or_default(),
            "meets_min_coverage",
        )?;
        let pass = parse_bool(
            record.get(schema.per_sample_qc_pass).unwrap_or_default(),
            "per_sample_qc_pass",
        )?;
        let normalized = cov
            .filter(|x| x.is_finite())
            .map(|x| x.ln_1p() - baselines[index]);
        current.push(Call {
            record,
            sample_index: index,
            cov,
            value,
            min_cov,
            pass,
            normalized_log_cov: normalized,
        });
        current_site = Some(site);
    }
    if let Some(site) = current_site {
        if process_site(
            &current,
            &site,
            &schema,
            args.total_samples,
            required,
            args.min_methylation_mad,
            &mut outputs,
            &mut scratch,
        )? {
            kept += 1;
        }
        processed += 1;
    }
    for writer in [
        &mut outputs.long,
        &mut outputs.qc,
        &mut outputs.metadata,
        &mut outputs.raw,
        &mut outputs.int,
    ] {
        writer.flush()?;
    }
    if temp_dir.exists() {
        fs::remove_dir_all(&temp_dir)?;
    }
    eprintln!(
        "Kept {kept} / {processed} sites after cohort-level QC for {}",
        args.chromosome
    );
    Ok(())
}
