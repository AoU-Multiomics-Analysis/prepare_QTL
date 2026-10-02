import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from collections import Counter
from reference_filter import filter_reference

def keep_cell(row):
    if row['matrix_header']=='CD8_T':
        return row['Monaco_fine_comparison']=='compatible'
    if row['matrix_header']=='Plasma':
        return not(row['Monaco_main_comparison']==row['Monaco_fine_comparison']=='uncertain')
    return True

def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    for key in ['reference','selected','predictions','gene-map','out-dir']:
        p.add_argument('--'+key,required=True)
    a=p.parse_args()
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    for path in [a.reference,a.selected,a.predictions,a.gene_map]:
        if '://' in path or not Path(path).is_file():raise ValueError('Unreadable local input')
    selected=read(a.selected);predictions=read(a.predictions)
    ids=[r['cell_id'] for r in selected]
    assert len(ids)==len(set(ids))==65191 and ids==[r['cell_id'] for r in predictions]
    assert all(x['matrix_header']==y['matrix_header'] for x,y in zip(selected,predictions))
    assert sha(a.reference)=='954e5e1e343ff241075c73dc67a18a81e3b9936830593a7f77d0a8ba5ead5c5b'
    mask=[keep_cell(r) for r in predictions]
    retained=[r for r,k in zip(selected,mask) if k]
    counts=dict(Counter(r['matrix_header'] for r in retained))
    assert len(retained)==63445 and counts['CD8_T']==262 and counts['Plasma']==522
    print('Filter reference: CD8 requires a retained fine CD8 call; plasma excludes calls uncertain in both. Other labels unchanged.',flush=True)
    validation=filter_reference(a.reference,out/'filtered_reference.tsv.gz',mask,[r['matrix_header'] for r in selected])
    gene_map=read(a.gene_map)
    sex_symbols={r['regulator_name_gencode'] for r in gene_map if r['regulator_chromosome'] in {'chrX','chrY','X','Y'}}
    genes=[]
    with gzip.open(out/'filtered_reference.tsv.gz','rt') as f:
        header=f.readline().rstrip('\n').split('\t')
        assert header[1:]==[r['matrix_header'] for r in retained]
        for line in f:genes.append(line.split('\t',1)[0])
    assert not(set(genes)&sex_symbols)
    assert not(set(genes)&{'HBB','HBA1','HBA2','KDM5D','XIST'})
    with (out/'input_gene_symbols.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['gene_symbol']);w.writerows([[g] for g in genes])
    for name,data in [('selected_cells.tsv',retained),('excluded_cells.tsv',[r for r,k in zip(predictions,mask) if not k])]:
        with (out/name).open('w',newline='') as f:
            if name=='selected_cells.tsv':
                fields=['parent_reference_matrix_column']+[k for k in data[0] if k!='matrix_column']+['matrix_column']
                data=[dict({k:v for k,v in r.items() if k!='matrix_column'},parent_reference_matrix_column=r['matrix_column'],matrix_column=i+2) for i,r in enumerate(data)]
            else:fields=list(data[0])
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    donors=Counter((r['matrix_header'],r['donor_id']) for r in retained)
    with (out/'counts_by_donor.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['cell_type','donor','retained_cells']);w.writerows([list(k)+[v] for k,v in sorted(donors.items())])
    validation.update(source=str(Path(a.reference).resolve()),source_sha256=sha(a.reference),gzip_sha256=sha(out/'filtered_reference.tsv.gz'),
        cells_before=65191,cells_excluded=1746,cell_type_counts=counts,sex_chromosome_genes_found=[],globin_exclusions_preserved=True,
        labels_reassigned=False,other_groups_filtered=False,CD8_rule='retained Monaco fine CD8 call',Plasma_rule='exclude both-model uncertainty')
    (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation,indent=2),flush=True)

if __name__=='__main__':main()
