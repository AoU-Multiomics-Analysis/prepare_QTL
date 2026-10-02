import gzip, hashlib, os
from pathlib import Path

def filter_reference(source, destination, keep, labels):
    keep=list(keep); labels=list(labels)
    if len(keep)!=len(labels) or not any(keep):
        raise ValueError('Invalid retention mask')
    dest=Path(destination); partial=dest.with_suffix(dest.suffix+'.partial')
    retained=[i+1 for i,k in enumerate(keep) if k]
    result={'cells':len(retained),'genes':0,'zero_genes_removed':[]}
    expected=hashlib.sha256()
    def write_line(out,line):
        payload=line.encode(); out.write(payload); expected.update(payload)
    with gzip.open(source,'rt') as src, gzip.open(partial,'wb',compresslevel=1) as out:
        header=src.readline().rstrip('\n').split('\t')
        if header[1:]!=labels: raise ValueError('Cell labels or order do not match the manifest')
        write_line(out,'\t'.join([header[0]]+[header[i] for i in retained])+'\n')
        for line in src:
            values=line.rstrip('\n').split('\t')
            if len(values)!=len(header): raise ValueError('Invalid reference row width')
            chosen=[values[i] for i in retained]
            if any(not v.isdecimal() for v in chosen): raise ValueError('Counts must be nonnegative integers')
            if not any(int(v)>0 for v in chosen):
                result['zero_genes_removed'].append(values[0]); continue
            write_line(out,'\t'.join([values[0]]+chosen)+'\n'); result['genes']+=1
    check=hashlib.sha256()
    with gzip.open(partial,'rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''): check.update(block)
    if check.hexdigest()!=expected.hexdigest(): raise ValueError('Output readback failed')
    os.replace(partial,dest)
    result['uncompressed_sha256']=expected.hexdigest()
    return result
