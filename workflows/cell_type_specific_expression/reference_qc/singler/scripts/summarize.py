import argparse
import collections
import csv
import gzip
import json
from pathlib import Path
from comparison import compare

def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def write(path, rows, fields=None):
    opener = gzip.open if str(path).endswith('.gz') else open
    if fields is None:
        fields = list(rows[0])
    with opener(path, 'wt', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def main():
    p = argparse.ArgumentParser()
    for key in ['predictions', 'flagged', 'clc-review', 'markers', 'out-dir']:
        p.add_argument('--' + key, required=True)
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = read(args.predictions)
    flagged_rows = read(args.flagged)
    flagged = {x['cell_id'] for x in flagged_rows if x['CD8_top10pct_S100'] == 'True'}
    clc_rows = read(args.clc_review)
    clc = {x['cell_id'] for x in clc_rows}
    markers = {x['cell_id']: x for x in read(args.markers)}
    assert len(rows) == 65191 and len({x['cell_id'] for x in rows}) == 65191
    assert set(markers) == {x['cell_id'] for x in rows}
    assert len(flagged) == 226 and len(clc) == 6
    models = ['Monaco_main', 'Monaco_fine', 'HPCA_main']
    counters = collections.Counter()
    for row in rows:
        row['flagged_S100_226'] = row['cell_id'] in flagged
        row['CLC_high_six'] = row['cell_id'] in clc
        statuses = []
        for model in models:
            label = row[model + '_label']
            pruned = row[model + '_pruned_label']
            if pruned == 'NA':
                pruned = ''
            status = compare(row['matrix_header'], label, pruned, model)
            row[model + '_comparison'] = status
            statuses.append(status)
            counters[(model, row['matrix_header'], row['donor_id'], label, status)] += 1
        row['review_disagreement'] = 'disagreement' in statuses
        row['review_uncertain'] = 'uncertain' in statuses
    fields = list(rows[0])
    write(out / 'all_cell_label_review.tsv.gz', rows, fields)
    write(out / 'label_disagreements.tsv.gz', [x for x in rows if x['review_disagreement']], fields)
    write(out / 'uncertain_calls.tsv.gz', [x for x in rows if x['review_uncertain']], fields)
    focused = [dict(x, **{k:v for k,v in markers[x['cell_id']].items() if k != 'cell_id'}) for x in rows if x['CLC_high_six'] or x['flagged_S100_226']]
    write(out / 'focused_CD8_review.tsv', focused)
    write(out / 'CLC_high_six_review.tsv', [x for x in focused if x['CLC_high_six']])
    write(out / 'retained_S100_51_review.tsv', [x for x in focused if x['flagged_S100_226']])
    stats = [dict(model=k[0], current_label=k[1], donor=k[2], predicted_label=k[3], comparison=k[4], cells=v)
             for k,v in sorted(counters.items())]
    write(out / 'label_counts_by_donor.tsv', stats)
    totals = collections.Counter()
    for k,v in counters.items():
        totals[(k[0],k[1],k[4])] += v
    write(out / 'comparison_totals.tsv', [dict(model=k[0], current_label=k[1], comparison=k[2], cells=v) for k,v in sorted(totals.items())])
    retained_by_id = {x['cell_id']: x for x in rows}
    all_flagged = []
    for row in flagged_rows:
        if row['cell_id'] not in flagged:
            continue
        if row['cell_id'] in retained_by_id:
            entry = dict(retained_by_id[row['cell_id']])
            entry['SingleR_review_status'] = 'reviewed_retained_singlet'
        else:
            entry = {k:'' for k in fields}
            for k in ['cell_id','donor_id','matrix_header','10X_run']:
                entry[k] = row.get(k,'')
            entry['SingleR_review_status'] = 'not_tested_removed_scDblFinder_doublet'
        all_flagged.append(entry)
    write(out / 'original_S100_226_review_status.tsv', all_flagged)
    summary = dict(reviewed_cells=len(rows), captures=len({x['10X_run'] for x in rows}), donors=len({x['donor_id'] for x in rows}),
        cells_with_any_disagreement=sum(x['review_disagreement'] for x in rows),
        cells_with_any_uncertain_call=sum(x['review_uncertain'] for x in rows),
        flagged_S100_retained=sum(x['flagged_S100_226'] for x in rows), CLC_high_reviewed=sum(x['CLC_high_six'] for x in rows),
        current_labels_changed=False, reference_rebuilt=False)
    assert summary['flagged_S100_retained'] == 51 and summary['CLC_high_reviewed'] == 6
    (out / 'validation.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary), flush=True)

if __name__ == '__main__':
    main()
