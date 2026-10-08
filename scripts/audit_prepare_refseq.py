#!/usr/bin/env python3
"""Offline RefSeq audit. No name guessing, accession-version fallback or species downsampling."""
import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from emito import cli


def op(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def fasta(path):
    name = None; parts = []
    with op(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: yield name, ''.join(parts).upper()
                name = line[1:].split()[0]; parts = []
            else:
                if name is None: raise ValueError('Sequence before header')
                parts.append(''.join(line.split()))
    if name is not None: yield name, ''.join(parts).upper()


def ancestor(tid, rank, nodes):
    seen = set()
    while tid not in seen and tid in nodes:
        seen.add(tid)
        parent, r = nodes[tid]
        if r == rank: return tid
        if parent == tid: break
        tid = parent
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--mapping', type=Path, help='NCBI accession2taxid(.gz), or accession_species.tsv; exact accession.version matches only')
    p.add_argument('--nodes-dmp', type=Path)
    p.add_argument('--names-dmp', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--prepare', action='store_true', help='Additionally split ALL QC-passing records and write metadata')
    p.add_argument('--min-genome-length', type=int, default=10000)
    a = p.parse_args()
    if a.prepare and not (a.mapping and a.nodes_dmp and a.names_dmp):
        p.error('--prepare requires --mapping, --nodes-dmp and --names-dmp')
    if a.mapping and not (a.nodes_dmp and a.names_dmp): p.error('Taxonomy files required with mapping')
    if a.min_genome_length < 1: p.error('Minimum length must be positive')
    if a.output.exists(): p.error('Output exists; specify a NEW directory')
    info = {}; counts = Counter()
    for acc, seq in fasta(a.fasta):
        if acc in info: raise ValueError(f'Duplicate accession in input: {acc}; not silently combined')
        if not seq or not seq.isascii() or not seq.isalpha(): raise ValueError(f'Invalid ungapped genome: {acc}')
        info[acc] = dict(accession=acc, length=len(seq), non_atcg_fraction=cli.fraction_bad(seq))
        counts['fasta_records'] += 1
        counts['records_below_length_threshold'] += len(seq) < a.min_genome_length
    if not info: raise ValueError('Empty FASTA')
    mapping = {}; nodes = {}; names = {}
    if a.mapping:
        with op(a.mapping) as f:
            reader = csv.DictReader(f, delimiter='\t')
            fields = reader.fieldnames or []
            acc_key = next((k for k in ['accession.version','accession_id','accession'] if k in fields), None)
            tid_key = next((k for k in ['original_taxid','taxid'] if k in fields), None)
            if not acc_key or not tid_key: raise ValueError('Unsupported mapping columns: '+str(fields))
            for row in reader:
                acc = row[acc_key].strip()
                if acc not in info: continue
                tid = row[tid_key].strip()
                if not tid.isdigit() or int(tid) <= 0: continue
                if acc in mapping and mapping[acc] != int(tid): raise ValueError('Conflicting TaxIDs: '+acc)
                mapping[acc] = int(tid)
        nodes = cli.core.load_nodes(a.nodes_dmp)
        wanted = set()
        for acc, row in info.items():
            tid = mapping.get(acc)
            sp = ancestor(tid, 'species', nodes)
            ge = ancestor(tid, 'genus', nodes)
            fa = ancestor(tid, 'family', nodes)
            row.update(original_taxid=tid, species_taxid=sp, genus_taxid=ge, family_taxid=fa)
            if sp: wanted.add(sp)
        names = cli.core.load_names(a.names_dmp, wanted)
    groups = defaultdict(list); retained = []
    for acc, row in info.items():
        sp = row.get('species_taxid')
        if sp: groups[sp].append(acc)
        row['species_name'] = names.get(sp, '')
        reasons = []
        if not a.mapping: reasons.append('taxonomy_not_checked')
        elif acc not in mapping: reasons.append('exact_accession_mapping_missing')
        elif not sp: reasons.append('species_ancestor_missing')
        elif not row.get('genus_taxid'): reasons.append('genus_ancestor_missing')
        if row['length'] < a.min_genome_length: reasons.append('below_length_threshold')
        row['status'] = ';'.join(reasons) if reasons else 'eligible'
        if not reasons: retained.append(acc)
    counts['mapped_records'] = len(mapping)
    counts['species_resolved'] = len(groups)
    counts['species_with_one_record'] = sum(len(v)==1 for v in groups.values())
    counts['species_with_multiple_records'] = sum(len(v)>1 for v in groups.values())
    counts['records_without_species_assignment'] = sum(not r.get('species_taxid') for r in info.values())
    counts['eligible_records'] = len(retained)
    counts['eligible_species'] = len({info[x]['species_taxid'] for x in retained})
    counts['minimum_length'] = min(r['length'] for r in info.values())
    counts['maximum_length'] = max(r['length'] for r in info.values())
    a.output.mkdir(parents=True)
    cli.write_tsv(a.output/'record_audit.tsv', ['accession','length','non_atcg_fraction','original_taxid',
        'species_taxid','species_name','genus_taxid','family_taxid','status'],info.values())
    cli.write_tsv(a.output/'species_counts.tsv',['species_taxid','species_name','records','accessions'],
        [dict(species_taxid=s,species_name=names.get(s,''),records=len(v),accessions=','.join(v))
         for s,v in sorted(groups.items(), key=lambda x:(-len(x[1]),x[0]))])
    cli.write_tsv(a.output/'summary.tsv',['metric','value'],[dict(metric=k,value=v) for k,v in counts.items()])
    if a.prepare:
        selected = set(retained)
        for acc, seq in fasta(a.fasta):
            if acc in selected:
                cli.write_genome_fasta(a.output/'prepared_input/fasta'/f'{acc}.fasta',[(acc,seq)])
        cli.write_tsv(a.output/'prepared_input/metadata.tsv',['accession_id','species_name','taxid','subgroup_label'],
            [dict(accession_id=acc,species_name=info[acc]['species_name'],taxid=mapping[acc],subgroup_label='') for acc in retained])
    cli.json_write(a.output/'audit_config.json',vars(a))
    (a.output/'COMPLETE').write_text('audit complete\n')
    print(json.dumps(dict(counts),indent=2))
    if not a.mapping: print('Species multiplicity UNKNOWN: provide an accession-to-TaxID mapping and taxonomy.')
    elif counts['records_without_species_assignment']: print('WARNING: unresolved records; species counts cover resolved records only.')
    print('Reports:', a.output)

if __name__ == '__main__':
    main()
