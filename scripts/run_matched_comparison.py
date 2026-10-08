#!/usr/bin/env python3
"""Species-only matched comparison: access -> two branches, without/with collapse."""
import argparse
import csv
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from emito import cli


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--metadata',type=Path)
    p.add_argument('--fasta-dir',type=Path)
    p.add_argument('--existing-standard',type=Path,help='Reuse COMPLETE fixed-version prepare/taxa/access from standard panel; no recomputation')
    p.add_argument('--names-dmp',type=Path)
    p.add_argument('--nodes-dmp',type=Path)
    p.add_argument('--mafft',default='mafft')
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--no-align',action='store_true',help='Skip MAFFT AND strand/origin normalization')
    p.add_argument('--min-genome-length',type=int,default=None,help='Minimum length in bp for fresh inputs (default: 10000)')
    p.add_argument('--require-single-reference',action='store_true',help='Fail unless metadata has exactly one record per species and no subgroup labels')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.existing_standard and a.min_genome_length is not None:
        p.error('Cannot change the length filter when reusing a standard panel; prepare a new panel first')
    if a.min_genome_length is None: a.min_genome_length = 10000
    if a.output.exists(): p.error('Use a NEW output directory')
    if a.existing_standard:
        for stage in ['01_prepare','02_taxa_generate','03_taxa_access']:
            if not (a.existing_standard/stage/'COMPLETE').is_file(): p.error('Incomplete standard stage: '+stage)
    elif not all([a.metadata,a.fasta_dir,a.names_dmp,a.nodes_dmp]): p.error('New inputs require metadata, fasta-dir, names-dmp, nodes-dmp')
    if a.min_genome_length < 1: p.error('Minimum genome length must be positive')
    if a.require_single_reference:
        if a.existing_standard: p.error('Single-reference check requires metadata inputs')
        from audit_prepare_refseq import ancestor
        nodes=cli.core.load_nodes(a.nodes_dmp)
        species=Counter(); accessions=set()
        for row in cli.read_tsv(a.metadata):
            if row.get('subgroup_label','').strip(): p.error('Single-reference metadata must have empty subgroup labels')
            acc=row['accession_id'].strip()
            if acc in accessions: p.error('Duplicate accession: '+acc)
            accessions.add(acc)
            sp=ancestor(int(row['taxid']),'species',nodes)
            if sp is None: p.error('Missing species ancestor: '+acc)
            species[sp]+=1
        duplicates={s:n for s,n in species.items() if n!=1}
        if not species or duplicates: p.error('Expected one reference per species; duplicates: '+str(duplicates))
        print(f'Confirmed {len(accessions):,} selected references for {len(species):,} species',flush=True)
        del nodes
    a.output.mkdir(parents=True)
    cli.json_write(a.output/'comparison_config.json',vars(a))
    def run(command,*args):
        cmd=[sys.executable,str(Path(__file__).with_name('eMito-'+command+'.py'))]+list(map(str,args))
        print('RUN:', ' '.join(cmd),flush=True)
        subprocess.run(cmd,check=True)
    def subset(src,dest):
        targets=cli.input_targets([src]); rows=[]
        dest.mkdir()
        for row in targets:
            if row['rank']=='species':
                r=dict(row); r['fasta']=str(row['path'].resolve())
                r['sequence_count']=sum(1 for _ in cli.core.iter_fasta(row['path']))
                rows.append(r)
        if not rows: raise ValueError('No species probe sets: '+str(src))
        cli.finish_bundle(dest,rows)
    prep=a.output/'01_prepare'
    taxa=a.output/'02_taxa_generate'
    access=a.output/'03_taxa_access'
    if a.existing_standard:
        prep=a.existing_standard/'01_prepare'
        subset(a.existing_standard/'02_taxa_generate',taxa)
        subset(a.existing_standard/'03_taxa_access',access)
    else:
        run('prepare','--metadata',a.metadata,'--fasta-dir',a.fasta_dir,'--names-dmp',a.names_dmp,
            '--nodes-dmp',a.nodes_dmp,'--mafft',a.mafft,'--threads',a.threads,
            '--min-genome-length',a.min_genome_length,
            '--no-align' if a.no_align else '--align','--output',prep)
        # Keep original routing, but disable subgroup output for the matched comparison.
        _,data=cli.prepared(prep)
        view=a.output/'01_species_view'; view.mkdir()
        for g in data['genomes']:
            g['group_eligible']=False
            for key in ['kmer_file','probe_file','normalized_fasta','organized_fasta']:
                g[key]=str((prep/g[key]).resolve())
        cli.json_write(view/'prepared.json',data)
        (view/'COMPLETE').write_text('Species-only view; prepared sequences unchanged\n')
        run('taxa-generate','--prepared',view,'--output',taxa)
        run('access','--inputs',taxa,'--output',access)
    run('merge','--inputs',access,'--output',a.output/'04_no_collapse_merge')
    run('collapse','--inputs',access,'--output',a.output/'05_collapse')
    run('merge','--inputs',a.output/'05_collapse','--output',a.output/'06_collapse_merge')
    _,data=cli.prepared(prep)
    members=[g for g in data['genomes'] if g['taxa_eligible']]
    rows=[dict(metric='input_genomes_species_route',value=len(members))]
    for key in ['species_taxid','genus_taxid']:
        rows.append(dict(metric=key+'_count',value=len({g[key] for g in members})))
    families={g['lineage'].get('family') for g in members}; families.discard(None)
    rows.append(dict(metric='assignable_families',value=len(families)))
    rows.append(dict(metric='genomes_without_family',value=sum('family' not in g['lineage'] for g in members)))
    for stage in [taxa,access,a.output/'05_collapse']:
        targets=cli.input_targets([stage])
        rows.append(dict(metric=stage.name+'_targets',value=len(targets)))
        rows.append(dict(metric=stage.name+'_records',value=sum(sum(1 for _ in cli.core.iter_fasta(t['path'])) for t in targets)))
    for stage in ['04_no_collapse_merge','06_collapse_merge']:
        summary=json.loads((a.output/stage/'merge_summary.json').read_text())
        rows.extend(dict(metric=stage+'_'+k,value=v) for k,v in summary.items())
    cli.write_tsv(a.output/'comparison_summary.tsv',['metric','value'],rows)
    (a.output/'COMPLETE').write_text('Both comparison branches completed\n')
    print('SUCCESS:',a.output,flush=True)

if __name__=='__main__':
    main()
