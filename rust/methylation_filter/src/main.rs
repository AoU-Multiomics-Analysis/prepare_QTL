use clap::Parser;
use csv::{ReaderBuilder, StringRecord, WriterBuilder};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{sync_channel, SyncSender};
use std::thread::{self, JoinHandle};
use thiserror::Error;

const REQUIRED_COLUMNS: [&str; 6] = ["#chrom", "begin", "end", "mod_score", "type", "cov"];
const OUTPUT_BATCH_ROWS: usize = 8192;
const OUTPUT_QUEUE_CAPACITY: usize = 4;

#[derive(Parser, Debug)]
#[command(about = "Stream pb-CpG methylation calls through per-sample QC and autosome splitting")]
struct Args {
    #[arg(long)]
    input_manifest: PathBuf,
    #[arg(long)]
    output_prefix: String,
    #[arg(long, default_value_t = 10.0)]
    min_coverage: f64,
    #[arg(long, default_value = "X|Y|M|_")]
    filter_chroms: String,
    #[arg(long, default_value_t = 3.0)]
    fence_k: f64,
    #[arg(long, default_value = "chr")]
    autosome_prefix: String,
    #[arg(long, default_value_t = 1)]
    num_threads: usize,
}

#[derive(Error, Debug)]
enum FilterError {
    #[error("{0}")]
    Message(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Csv(#[from] csv::Error),
}

type Result<T> = std::result::Result<T, FilterError>;

#[derive(Clone)]
struct InputColumns {
    chrom: usize,
    begin: usize,
    end: usize,
    mod_score: usize,
    call_type: usize,
    coverage: usize,
}

struct SampleStats {
    input_rows: u64,
    rows_after_chrom_filter: u64,
    call_type: String,
    median_coverage: f64,
    extreme_cutoff: f64,
}

fn open_reader(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).map_err(|error| {
        FilterError::Message(format!("Cannot open {}: {error}", path.display()))
    })?;
    if path
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
    {
        Ok(Box::new(MultiGzDecoder::new(BufReader::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn reader_for_bed(path: &Path) -> Result<(csv::Reader<Box<dyn Read>>, InputColumns)> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_reader(open_reader(path)?);
    let mut header = StringRecord::new();
    loop {
        if !reader.read_record(&mut header)? {
            return Err(FilterError::Message(format!(
                "Missing #chrom header in {}",
                path.display()
            )));
        }
        if header.get(0) == Some("#chrom") {
            break;
        }
    }
    let indices = REQUIRED_COLUMNS.map(|column| {
        header
            .iter()
            .position(|candidate| candidate == column)
            .ok_or_else(|| {
                FilterError::Message(format!(
                    "Missing required column '{column}' in {}",
                    path.display()
                ))
            })
    });
    let [chrom, begin, end, mod_score, call_type, coverage] = indices;
    Ok((
        reader,
        InputColumns {
            chrom: chrom?,
            begin: begin?,
            end: end?,
            mod_score: mod_score?,
            call_type: call_type?,
            coverage: coverage?,
        },
    ))
}

fn parse_coverage(value: &str, path: &Path) -> Result<Option<f64>> {
    if value.is_empty() || value == "NA" || value == "." {
        return Ok(None);
    }
    value.parse::<f64>().map(Some).map_err(|_| {
        FilterError::Message(format!("Column 'cov' is not numeric in {}", path.display()))
    })
}

fn type7_quantile(sorted: &[f64], probability: f64) -> f64 {
    if sorted.len() == 1 {
        return sorted[0];
    }
    let index = (sorted.len() - 1) as f64 * probability;
    let lower = index.floor() as usize;
    let upper = index.ceil() as usize;
    sorted[lower] + (index - lower as f64) * (sorted[upper] - sorted[lower])
}

fn scan_sample(path: &Path, filter_chroms: Option<&Regex>, fence_k: f64) -> Result<SampleStats> {
    let (mut reader, columns) = reader_for_bed(path)?;
    let mut record = StringRecord::new();
    let mut input_rows = 0_u64;
    let mut retained_rows = 0_u64;
    let mut coverages = Vec::<f64>::new();
    let mut log_coverages = Vec::<f64>::new();
    let mut call_types = HashSet::<String>::new();
    let mut sites = HashSet::<(u32, u64, u64)>::new();
    let mut chromosome_ids = HashMap::<String, u32>::new();
    while reader.read_record(&mut record)? {
        input_rows += 1;
        let chrom = record.get(columns.chrom).unwrap_or_default();
        if filter_chroms.is_some_and(|pattern| pattern.is_match(chrom)) {
            continue;
        }
        retained_rows += 1;
        let begin = record
            .get(columns.begin)
            .unwrap_or_default()
            .parse::<u64>()
            .map_err(|_| {
                FilterError::Message(format!(
                    "Column 'begin' is not an unsigned integer in {}",
                    path.display()
                ))
            })?;
        let end = record
            .get(columns.end)
            .unwrap_or_default()
            .parse::<u64>()
            .map_err(|_| {
                FilterError::Message(format!(
                    "Column 'end' is not an unsigned integer in {}",
                    path.display()
                ))
            })?;
        let next_chromosome_id = chromosome_ids.len() as u32;
        let chromosome_id = *chromosome_ids
            .entry(chrom.to_owned())
            .or_insert(next_chromosome_id);
        if !sites.insert((chromosome_id, begin, end)) {
            return Err(FilterError::Message(format!(
                "Found duplicated #chrom/begin/end site(s) in {}. Aggregate duplicate calls before merging so a site is counted once per sample.",
                path.display()
            )));
        }
        let call_type = record.get(columns.call_type).unwrap_or_default();
        call_types.insert(call_type.to_owned());
        if let Some(coverage) =
            parse_coverage(record.get(columns.coverage).unwrap_or_default(), path)?
        {
            if coverage.is_finite() {
                coverages.push(coverage);
            }
            if coverage.is_finite() && coverage > 0.0 {
                log_coverages.push(coverage.log10());
            }
        }
    }
    if retained_rows == 0 {
        return Err(FilterError::Message(format!(
            "No rows remain after chromosome filtering in {}",
            path.display()
        )));
    }
    if call_types.len() != 1 {
        return Err(FilterError::Message(format!(
            "Expected one pb-CpG-tools 'type' per input file, but found: {}. Use one .combined.bed.gz file per sample.",
            call_types.into_iter().collect::<Vec<_>>().join(", ")
        )));
    }
    coverages.sort_by(|left, right| left.total_cmp(right));
    let median_coverage = if coverages.is_empty() {
        f64::NAN
    } else {
        type7_quantile(&coverages, 0.5)
    };
    if !median_coverage.is_finite() || median_coverage <= 0.0 {
        return Err(FilterError::Message(format!(
            "Median coverage must be positive after chromosome filtering in {}",
            path.display()
        )));
    }
    log_coverages.sort_by(|left, right| left.total_cmp(right));
    let extreme_cutoff = if log_coverages.is_empty() {
        f64::INFINITY
    } else {
        let q25 = type7_quantile(&log_coverages, 0.25);
        let q75 = type7_quantile(&log_coverages, 0.75);
        10_f64.powf(q75 + fence_k * (q75 - q25))
    };
    Ok(SampleStats {
        input_rows,
        rows_after_chrom_filter: retained_rows,
        call_type: call_types.into_iter().next().unwrap(),
        median_coverage,
        extreme_cutoff,
    })
}

fn format_number(value: f64) -> String {
    if value.is_nan() {
        "NA".to_owned()
    } else if value.is_infinite() {
        if value.is_sign_positive() {
            "Inf".to_owned()
        } else {
            "-Inf".to_owned()
        }
    } else {
        value.to_string()
    }
}

fn output_path(prefix: &str, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{prefix}{suffix}"))
}

struct CompressionWorker {
    sender: SyncSender<Option<(usize, Vec<u8>)>>,
    join: JoinHandle<Result<()>>,
}

struct ParallelWriters {
    buffers: Vec<csv::Writer<Vec<u8>>>,
    buffered_rows: Vec<usize>,
    workers: Vec<CompressionWorker>,
}

fn new_batch_writer() -> csv::Writer<Vec<u8>> {
    WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(Vec::with_capacity(1024 * 1024))
}

impl ParallelWriters {
    fn new(output_prefix: &str, header: [&str; 12], task_cpus: usize) -> Result<Self> {
        // Reserve one task CPU for parsing and route each chromosome to one
        // deterministic compressor, preserving record order within its file.
        let worker_count = task_cpus.saturating_sub(1).clamp(1, 22);
        let mut workers = Vec::with_capacity(worker_count);
        for worker_index in 0..worker_count {
            let (sender, receiver) =
                sync_channel::<Option<(usize, Vec<u8>)>>(OUTPUT_QUEUE_CAPACITY);
            let prefix = output_prefix.to_owned();
            let join = thread::spawn(move || -> Result<()> {
                let mut outputs: Vec<Option<BufWriter<GzEncoder<File>>>> =
                    (0..22).map(|_| None).collect();
                for chromosome_index in 0..22 {
                    if chromosome_index % worker_count == worker_index {
                        let path = output_path(
                            &prefix,
                            &format!(
                                ".methylation.autosome{:02}.per_sample_qc.long.tsv.gz",
                                chromosome_index + 1
                            ),
                        );
                        let encoder = GzEncoder::new(File::create(path)?, Compression::new(3));
                        outputs[chromosome_index] =
                            Some(BufWriter::with_capacity(1024 * 1024, encoder));
                    }
                }
                while let Some((chromosome_index, bytes)) = receiver.recv().map_err(|_| {
                    FilterError::Message(
                        "Parallel compression worker stopped unexpectedly".to_owned(),
                    )
                })? {
                    let output = outputs
                        .get_mut(chromosome_index)
                        .and_then(|output| output.as_mut())
                        .ok_or_else(|| {
                            FilterError::Message("Invalid chromosome compression route".to_owned())
                        })?;
                    output.write_all(&bytes)?;
                }
                for output in outputs.into_iter().flatten() {
                    output
                        .into_inner()
                        .map_err(|error| error.into_error())?
                        .finish()?;
                }
                Ok(())
            });
            workers.push(CompressionWorker { sender, join });
        }

        let mut buffers = Vec::with_capacity(22);
        for _ in 0..22 {
            let mut buffer = new_batch_writer();
            buffer.write_record(header)?;
            buffers.push(buffer);
        }
        Ok(Self {
            buffers,
            buffered_rows: vec![0; 22],
            workers,
        })
    }

    fn flush_batch(&mut self, chromosome_index: usize) -> Result<()> {
        let buffer = std::mem::replace(&mut self.buffers[chromosome_index], new_batch_writer());
        let bytes = buffer.into_inner().map_err(|error| {
            FilterError::Message(format!(
                "Cannot finalize chromosome output batch: {}",
                error.error()
            ))
        })?;
        self.buffered_rows[chromosome_index] = 0;
        if !bytes.is_empty() {
            self.workers[chromosome_index % self.workers.len()]
                .sender
                .send(Some((chromosome_index, bytes)))
                .map_err(|_| {
                    FilterError::Message(
                        "Parallel compression worker stopped unexpectedly".to_owned(),
                    )
                })?;
        }
        Ok(())
    }

    fn write_record(&mut self, chromosome_index: usize, record: [&str; 12]) -> Result<()> {
        self.buffers[chromosome_index].write_record(record)?;
        self.buffered_rows[chromosome_index] += 1;
        if self.buffered_rows[chromosome_index] >= OUTPUT_BATCH_ROWS {
            self.flush_batch(chromosome_index)?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        for chromosome_index in 0..self.buffers.len() {
            self.flush_batch(chromosome_index)?;
        }
        for worker in &self.workers {
            worker.sender.send(None).map_err(|_| {
                FilterError::Message("Parallel compression worker stopped unexpectedly".to_owned())
            })?;
        }
        for worker in self.workers {
            worker.join.join().map_err(|_| {
                FilterError::Message("Parallel compression worker panicked".to_owned())
            })??;
        }
        Ok(())
    }
}

fn write_sample(
    sample_id: &str,
    file_path: &Path,
    stats: &SampleStats,
    filter_chroms: Option<&Regex>,
    min_coverage: f64,
    autosome_indices: &HashMap<String, usize>,
    writers: &mut ParallelWriters,
) -> Result<(u64, u64, u64, u64)> {
    let (mut reader, columns) = reader_for_bed(file_path)?;
    let mut record = StringRecord::new();
    let mut below_minimum = 0_u64;
    let mut extreme_total = 0_u64;
    let mut extreme_after_minimum = 0_u64;
    let mut passing = 0_u64;
    while reader.read_record(&mut record)? {
        let chrom = record.get(columns.chrom).unwrap_or_default();
        if filter_chroms.is_some_and(|pattern| pattern.is_match(chrom)) {
            continue;
        }
        let coverage = parse_coverage(record.get(columns.coverage).unwrap_or_default(), file_path)?;
        let coverage_pass = coverage.is_some_and(|value| value >= min_coverage);
        let extreme_pass = !coverage.is_some_and(|value| value >= stats.extreme_cutoff);
        if !coverage_pass {
            below_minimum += 1;
        }
        if !extreme_pass {
            extreme_total += 1;
        }
        if coverage_pass && !extreme_pass {
            extreme_after_minimum += 1;
        }
        let per_sample_pass = coverage_pass && extreme_pass;
        if per_sample_pass {
            passing += 1;
        }
        if let Some(&chromosome_index) = autosome_indices.get(chrom) {
            let cov_text = coverage
                .map(format_number)
                .unwrap_or_else(|| "NA".to_owned());
            let implied_cn = coverage
                .map(|value| format_number(2.0 * value / stats.median_coverage))
                .unwrap_or_else(|| "NA".to_owned());
            let extreme_flag = if extreme_pass { "ok" } else { "extreme" };
            let begin = record.get(columns.begin).unwrap_or_default();
            let end = record.get(columns.end).unwrap_or_default();
            let site_key = format!("{chrom}*{begin}*{end}");
            let meets_min_coverage = if coverage_pass { "TRUE" } else { "FALSE" };
            let per_sample_qc_pass = if per_sample_pass { "TRUE" } else { "FALSE" };
            writers.write_record(
                chromosome_index,
                [
                    chrom,
                    begin,
                    end,
                    record.get(columns.mod_score).unwrap_or_default(),
                    record.get(columns.call_type).unwrap_or_default(),
                    &cov_text,
                    &implied_cn,
                    extreme_flag,
                    &site_key,
                    sample_id,
                    meets_min_coverage,
                    per_sample_qc_pass,
                ],
            )?;
        }
    }
    Ok((below_minimum, extreme_total, extreme_after_minimum, passing))
}

fn main() -> Result<()> {
    let args = Args::parse();
    if !args.min_coverage.is_finite() || args.min_coverage < 0.0 {
        return Err(FilterError::Message(
            "--MinCoverage must be a non-negative number".to_owned(),
        ));
    }
    if !args.fence_k.is_finite() || args.fence_k < 0.0 {
        return Err(FilterError::Message(
            "--FenceK must be a non-negative number".to_owned(),
        ));
    }
    if args.num_threads < 1 {
        return Err(FilterError::Message(
            "--NumThreads must be at least 1".to_owned(),
        ));
    }
    eprintln!(
        "Streaming methylation filtering with {} task CPU(s)",
        args.num_threads
    );
    let filter_chroms = if args.filter_chroms.is_empty() {
        None
    } else {
        Some(Regex::new(&args.filter_chroms).map_err(|error| {
            FilterError::Message(format!("Invalid --FilterChroms regex: {error}"))
        })?)
    };
    let mut manifest_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&args.input_manifest)?;
    let manifest_headers = manifest_reader.headers()?.clone();
    let sample_index = manifest_headers
        .iter()
        .position(|column| column == "sample_id")
        .ok_or_else(|| {
            FilterError::Message("Input manifest must contain column 'sample_id'".to_owned())
        })?;
    let path_index = manifest_headers
        .iter()
        .position(|column| column == "file_path")
        .ok_or_else(|| {
            FilterError::Message("Input manifest must contain column 'file_path'".to_owned())
        })?;
    let manifest_dir = args
        .input_manifest
        .parent()
        .unwrap_or_else(|| Path::new("."));
    let autosomes: Vec<String> = (1..=22)
        .map(|index| format!("{}{}", args.autosome_prefix, index))
        .collect();
    let autosome_indices: HashMap<String, usize> = autosomes
        .iter()
        .enumerate()
        .map(|(index, chromosome)| (chromosome.clone(), index))
        .collect();
    let output_parent = Path::new(&args.output_prefix)
        .parent()
        .unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_parent)?;
    let header = [
        "#chrom",
        "begin",
        "end",
        "mod_score",
        "type",
        "cov",
        "implied_cn",
        "extreme_cov_flag",
        "site_key",
        "sample_id",
        "meets_min_coverage",
        "per_sample_qc_pass",
    ];
    let mut writers = ParallelWriters::new(&args.output_prefix, header, args.num_threads)?;
    let mut qc_rows = Vec::new();
    let mut sample_ids = HashSet::new();
    let mut reference_columns: Option<[String; 6]> = None;
    for (sample_number, row) in manifest_reader.records().enumerate() {
        let row = row?;
        let sample_id = row.get(sample_index).unwrap_or_default();
        let raw_path = row.get(path_index).unwrap_or_default();
        if sample_id.is_empty() || raw_path.is_empty() {
            return Err(FilterError::Message(
                "Input manifest contains an empty sample_id or file_path".to_owned(),
            ));
        }
        if !sample_ids.insert(sample_id.to_owned()) {
            return Err(FilterError::Message(
                "Each sample_id must occur exactly once in the input manifest".to_owned(),
            ));
        }
        let file_path = {
            let path = PathBuf::from(raw_path);
            if path.is_absolute() {
                path
            } else {
                manifest_dir.join(path)
            }
        };
        if !file_path.exists() {
            return Err(FilterError::Message(format!(
                "Input BED file does not exist: {}",
                file_path.display()
            )));
        }
        eprintln!(
            "[{}] Processing {}: {}",
            sample_number + 1,
            sample_id,
            file_path.display()
        );
        let stats = scan_sample(&file_path, filter_chroms.as_ref(), args.fence_k)?;
        let current_columns = REQUIRED_COLUMNS.map(str::to_owned);
        if let Some(reference) = &reference_columns {
            if reference != &current_columns {
                return Err(FilterError::Message(
                    "BED columns do not match the first input file".to_owned(),
                ));
            }
        } else {
            reference_columns = Some(current_columns);
        }
        let (below_minimum, extreme_total, extreme_after_minimum, passing) = write_sample(
            sample_id,
            &file_path,
            &stats,
            filter_chroms.as_ref(),
            args.min_coverage,
            &autosome_indices,
            &mut writers,
        )?;
        eprintln!(
            "  Input sites: {}; removed by chromosome filter: {}; evaluated for coverage: {}",
            stats.input_rows,
            stats.input_rows - stats.rows_after_chrom_filter,
            stats.rows_after_chrom_filter
        );
        eprintln!("  Per-sample thresholds: {below_minimum} fail MinCoverage; {extreme_after_minimum} fail extreme coverage after MinCoverage; {passing} pass both thresholds");
        qc_rows.push(vec![
            sample_id.to_owned(),
            raw_path.to_owned(),
            stats.input_rows.to_string(),
            stats.rows_after_chrom_filter.to_string(),
            (stats.input_rows - stats.rows_after_chrom_filter).to_string(),
            stats.call_type,
            format_number(stats.median_coverage),
            format_number(stats.extreme_cutoff),
            below_minimum.to_string(),
            extreme_total.to_string(),
            extreme_after_minimum.to_string(),
            passing.to_string(),
        ]);
    }
    if qc_rows.is_empty() {
        return Err(FilterError::Message(
            "Input manifest must contain at least one data row".to_owned(),
        ));
    }
    writers.finish()?;
    let qc_path = output_path(&args.output_prefix, ".methylation.sample_qc.tsv");
    let mut qc_writer = WriterBuilder::new().delimiter(b'\t').from_path(qc_path)?;
    qc_writer.write_record([
        "sample_id",
        "file_path",
        "n_input_rows",
        "n_rows_after_chrom_filter",
        "n_removed_by_chrom_filter",
        "pb_cpg_type",
        "median_cov",
        "extreme_coverage_cutoff",
        "n_below_min_coverage",
        "n_extreme_coverage",
        "n_extreme_coverage_after_min_coverage",
        "n_passing_per_sample_qc",
    ])?;
    for row in qc_rows {
        qc_writer.write_record(row)?;
    }
    qc_writer.flush()?;
    eprintln!("Wrote sample QC and 22 autosome-split QC-flagged call tables");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::type7_quantile;

    #[test]
    fn matches_r_type7_quantiles() {
        let values = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(type7_quantile(&values, 0.25), 1.75);
        assert_eq!(type7_quantile(&values, 0.5), 2.5);
        assert_eq!(type7_quantile(&values, 0.75), 3.25);
    }
}
