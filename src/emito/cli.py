"""Independent, composable eMito commands. No installed legacy eMito required."""
import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path

from . import __version__
from . import pipeline as core

RC = str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn', 'TGCAYRKMVHDBNtgcayrkmvhdbn')
core.reverse_complement = lambda s: s.translate(RC)[::-1]


def read_tsv(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def json_write(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2, default=str)+'\n', encoding='utf-8')


def new_output(args):
    p = args.output.resolve()
    if p.exists():
        raise core.PipelineError(f'Output already exists; choose a NEW directory: {p}')
    p.mkdir(parents=True)
    json_write(p/'config.json', dict(vars(args), version=__version__))
    return p


def fraction_bad(seq):
    return sum(b not in 'ACGT' for b in seq.upper()) / len(seq)


def reference_key(item):
    return (not item['accession'].startswith('NC_'), item['non_atcg_fraction'], item['accession'])


def circular_windows(seq, length, step):
    if len(seq) < length:
        return
    extended = seq + seq[:length-1]
    for start in range(0, len(seq), step):
        window = extended[start:start+length]
        if core.DNA_RE.fullmatch(window):
            yield start, window


def window_header(g, start, length, kind, reverse=False, alignment_start=None):
    """Normalized interval is unrolled (end may exceed genome length)."""
    n = g['length']
    end = start+length-1
    orig_start = core.original_coordinate(start, n, g['orientation'], g['rotation'])
    orig_end = core.original_coordinate(end, n, g['orientation'], g['rotation'])
    strand = 'ReverseComplement' if reverse else 'Forward'
    if reverse:
        s = (n-1-(end % n)) % n + 1
        e = (n-1-start) % n + 1
    else:
        s, e = orig_start, orig_end
    header = (f"{g['accession']}|start={s}|end={e}|strand={strand}"
              f"|normalized_start={start+1}|normalized_end={end+1}"
              f"|genome_length={n}|circular=true|wrap={int(end>=n)}"
              f"|coordinate_frame={g['coordinate_frame']}"
              f"|source_orientation={g['orientation']}|rotation={g['rotation']}|type={kind}")
    if alignment_start is not None:
        header += f'|alignment_start={alignment_start+1}'
    if reverse:
        header += ' [Reverse Complement]'
    return header


def write_genome_fasta(path, records):
    """Preserve ambiguity codes; ATCG filtering belongs to windows, not genomes."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='ascii') as out:
        for header, sequence in records:
            if not sequence or not re.fullmatch('[A-Z]+', sequence):
                raise core.PipelineError(f'Invalid ungapped genome: {header}')
            out.write(f'>{header}\n{sequence}\n')


def prepare_group(task):
    members, root, align, mafft, anchor, min_anchor, min_hits, length, pstep, kstep = task
    root = Path(root)
    reference = min(members, key=reference_key)
    group_dir = root/'alignments'/reference['alignment_group']
    group_dir.mkdir(parents=True)
    ref = core.read_single_fasta(Path(reference['organized_fasta']))[1]
    normalized = {}
    for g in members:
        seq = core.read_single_fasta(Path(g['organized_fasta']))[1]
        if not align or g['accession'] == reference['accession']:
            seqn, orientation, rotation, used_anchor, hits = seq, 'Forward', 0, 0, 0
        else:
            result = core.circular_normalize_to_reference(ref, seq, anchor, min_anchor, min_hits)
            seqn, orientation, rotation = result.sequence, result.orientation, result.rotation
            used_anchor, hits = result.anchor_length, result.anchor_hits
        g.update(orientation=orientation, rotation=rotation, reference_accession=reference['accession'],
                 anchor_length_used=used_anchor, anchor_hits=hits,
                 coordinate_frame=f"{reference['accession']}:{orientation}:{rotation}" if align else 'raw')
        g['normalized_fasta'] = f"normalized/{g['accession']}.fasta"
        write_genome_fasta(root/g['normalized_fasta'], [(g['accession'], seqn)])
        normalized[g['accession']] = seqn
    inp = group_dir/'normalized_input.fasta'
    write_genome_fasta(inp, sorted(normalized.items()))
    aln = group_dir/'alignment.fasta'
    if align and len(members)>1:
        with aln.open('w') as out, (group_dir/'mafft.log').open('w') as err:
            proc = subprocess.run([mafft, '--auto', '--thread', '1', str(inp)], stdout=out, stderr=err)
        if proc.returncode:
            raise core.PipelineError(f'MAFFT failed: {group_dir}')
        aligned = core.read_alignment(aln)
    else:
        aligned = normalized.copy()
        if align:
            shutil.copyfile(inp, aln)
        (group_dir/'mafft.log').write_text('Singleton: no MAFFT\n' if align else 'Alignment disabled\n')
    missing = sorted(set(normalized) - set(aligned))
    extra = sorted(set(aligned) - set(normalized))
    changed = sorted(a for a in set(normalized) & set(aligned)
                     if aligned[a].replace('-', '') != normalized[a])
    if missing or extra or changed:
        raise core.PipelineError(f'Alignment changed input sequences: {group_dir}; '
                                 f'missing={missing}; unexpected={extra}; changed={changed}')
    columns = core.aligned_reference_columns(aligned[reference['accession']]) if align else None
    for g in members:
        acc = g['accession']
        seq = normalized[acc]
        if align:
            _, offsets, _ = core.alignment_coordinate_maps(aligned[acc])
        def windows(step):
            if len(seq)<length:
                return
            if not align:
                for start, w in circular_windows(seq, length, step):
                    yield start, w, None
            else:
                for roffset in range(0, len(columns), step):
                    col = columns[roffset]
                    if aligned[acc][col] == '-':
                        continue
                    start = offsets[col]
                    w = (seq+seq[:length-1])[start:start+length]
                    if core.DNA_RE.fullmatch(w):
                        yield start,w,col
        g['probe_file'] = f'windows/{acc}.probe.fasta'
        g['kmer_file'] = f'windows/{acc}.kmer.merged.fasta'
        g['probe_count'] = core.write_fasta(root/g['probe_file'],
            ((window_header(g,s,length,'probe',alignment_start=c),w) for s,w,c in windows(pstep)))
        paths = [root/f'windows/{acc}.kmer.{suffix}.fasta' for suffix in ['forward','reverse_complement','merged']]
        g['kmer_forward_count'] = 0
        with paths[0].open('w') as f, paths[1].open('w') as r, paths[2].open('w') as m:
            for s,w,c in windows(kstep):
                ft = f">{window_header(g,s,length,'kmer',alignment_start=c)}\n{w}\n"
                rt = f">{window_header(g,s,length,'kmer',True,c)}\n{core.reverse_complement(w)}\n"
                f.write(ft); r.write(rt); m.write(ft+rt)
                g['kmer_forward_count'] += 1
    return members


def prepare(args):
    if args.circular_min_anchor_length > args.circular_anchor_length:
        raise core.PipelineError('Minimum anchor length exceeds initial length')
    mafft = core.resolve_executable(args.mafft,'MAFFT') if args.align else None
    nodes = core.load_nodes(args.nodes_dmp)
    genomes = core.build_genomes(core.load_metadata(args.metadata), args.fasta_dir, nodes,
                                  args.names_dmp, args.output.resolve()/'organized')
    genomes, qc = core.filter_genomes_by_length(genomes,args.min_genome_length)
    if not genomes:
        raise core.PipelineError('No genomes pass length QC')
    taxa, groups, _ = core.route_genomes_by_mode(genomes)
    ta = {g.accession for g in taxa}; ga = {g.accession for g in groups}
    bysp = defaultdict(list)
    for g in genomes: bysp[g.species_taxid].append(g)
    special = set()
    for sp, members in bysp.items():
        direct = [g for g in members if not g.subgroup_label]
        subgroup = [g for g in members if g.subgroup_label]
        if direct and len({g.subgroup_label for g in subgroup}) == 1 and all(g.taxid!=sp for g in subgroup):
            for g in subgroup:
                ga.add(g.accession); special.add(g.accession)
    root = new_output(args)
    core.write_genome_length_qc(root,qc)
    core.organize_genomes(genomes,'copy')
    grouped = defaultdict(list)
    for g in genomes:
        sequence = core.read_single_fasta(g.organized_fasta)[1]
        if not re.fullmatch('[A-Z]+',sequence):
            raise core.PipelineError(f'Use ungapped input FASTA (letters only): {g.accession}')
        groupid = core.make_group_target(g).target_id if g.subgroup_label else f'species_{g.species_taxid}'
        row = asdict(g)
        row.update(source_fasta=str(g.source_fasta.resolve()),organized_fasta=str(g.organized_fasta),
                   length=len(sequence),non_atcg_fraction=fraction_bad(sequence),alignment_group=groupid,
                   taxa_eligible=g.accession in ta, group_eligible=g.accession in ga,
                   extra_species_background=g.accession in special,
                   group_target=core.make_group_target(g).target_id if g.subgroup_label else '')
        grouped[groupid].append(row)
    tasks = [(grouped[k],str(root),args.align,mafft,args.circular_anchor_length,
              args.circular_min_anchor_length,args.circular_min_anchor_hits,args.window_length,
              args.probe_step,args.kmer_step) for k in sorted(grouped)]
    result=[]
    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        for i,members in enumerate(pool.map(prepare_group,tasks),1):
            result.extend(members)
            print(f'prepare: {i}/{len(tasks)} groups',flush=True)
    # Save complete ancestor maps so node selection needs no external taxdump.
    wanted=set()
    for row in result:
        lineage={}; tid=row['taxid']; visited=set()
        while tid not in visited:
            visited.add(tid)
            parent,rank=nodes[tid]
            if rank not in lineage: lineage[rank]=tid
            wanted.add(tid)
            if parent==tid: break
            tid=parent
        row['lineage']=lineage
        row['organized_fasta']=str(Path(row['organized_fasta']).relative_to(root))
    names=core.load_names(args.names_dmp,wanted)
    json_write(root/'prepared.json',dict(schema='emito.prepared.v1',version=__version__,
        window_length=args.window_length,probe_step=args.probe_step,kmer_step=args.kmer_step,
        aligned=args.align,names=names,genomes=result))
    write_tsv(root/'routing.tsv', ['accession','species_name','species_taxid','taxid','subgroup_label',
        'taxa_eligible','group_eligible','extra_species_background','reference_accession',
        'non_atcg_fraction','length','orientation','rotation','anchor_length_used','anchor_hits',
        'probe_count','kmer_forward_count'], result)
    (root/'COMPLETE').write_text('prepare\n')


def prepared(path):
    path=Path(path).resolve()
    if not (path/'COMPLETE').is_file(): raise core.PipelineError(f'Incomplete prepare output: {path}')
    data=json.loads((path/'prepared.json').read_text())
    if data['schema']!='emito.prepared.v1': raise core.PipelineError('Unknown prepare schema')
    return path,data


def finish_bundle(root,rows):
    write_tsv(root/'manifest.tsv',['target_id','rank','name','taxid','fasta','sequence_count'],rows)
    json_write(root/'bundle.json',dict(schema='emito.probes.v1',version=__version__))
    (root/'COMPLETE').write_text('probe bundle\n')


def bundle_row(target,rank,name,taxid,file,count):
    return dict(target_id=target,rank=rank,name=name,taxid=taxid,fasta=str(file),sequence_count=count)


def representative(root,members,args):
    return core.representative_kmers([root/g['kmer_file'] for g in members],args.representative_fraction,
                                    args.small_group_all_max)


def specificity_outputs(out, prep_root, members, allowed, key, rank, name, taxid):
    core.write_fasta(out/f'specific_kmers/{key}.fasta', ((h,s) for s,h in sorted(allowed.items())))
    files=[]
    for g in sorted(members,key=lambda x:x['accession']):
        file=out/f"intersections/{key}/{g['accession']}.fasta"
        core.intersect_probe_file(prep_root/g['probe_file'],set(allowed),file)
        files.append(file)
    file=Path('probe_sets')/f'{key}.fasta'
    count=core.merge_probe_files_deduplicated(files,out/file)
    return bundle_row(key,rank,name,taxid,file,count)


def taxa_generate(args):
    if not 0<args.representative_fraction<=1: raise core.PipelineError('Fraction must be in (0,1]')
    prep,data=prepared(args.prepared)
    out=new_output(args)
    species=defaultdict(list); groups=defaultdict(list)
    for g in data['genomes']:
        if g['taxa_eligible']: species[g['species_taxid']].append(g)
        if g['group_eligible']: groups[g['group_target']].append(g)
    # Owner sets are exact representative-set memberships, not raw-genome occurrence.
    sp_owners=defaultdict(set); genus_owners=defaultdict(set); threshold_rows=[]
    for sp,members in sorted(species.items()):
        rep,t=representative(prep,members,args)
        core.write_fasta(out/f'representative_kmers/species_{sp}.fasta',((h,s) for s,h in rep.items()))
        threshold_rows.append(dict(target_id=f'species_{sp}',genomes=len(members),required=t,kmers=len(rep)))
        for s in rep:
            sp_owners[s].add(sp); genus_owners[s].add(members[0]['genus_taxid'])
    rows=[]
    for sp,members in sorted(species.items()):
        rep=core.sequence_dictionary(out/f'representative_kmers/species_{sp}.fasta')
        onlysp={s:h for s,h in rep.items() if len(sp_owners[s])==1}
        onlygenus={s:h for s,h in rep.items() if len(genus_owners[s])==1}
        core.write_fasta(out/f'specific_kmers/species_{sp}.species_only.fasta',((h,s) for s,h in onlysp.items()))
        core.write_fasta(out/f'specific_kmers/species_{sp}.genus_restricted.fasta',((h,s) for s,h in onlygenus.items()))
        allowed=dict(onlygenus); allowed.update(onlysp)
        rows.append(specificity_outputs(out,prep,members,allowed,f'species_{sp}','species',members[0]['species_name'],sp))
    del genus_owners
    owners=defaultdict(set)
    for key,members in sorted(groups.items()):
        rep,t=representative(prep,members,args)
        core.write_fasta(out/f'representative_kmers/{key}.fasta',((h,s) for s,h in rep.items()))
        threshold_rows.append(dict(target_id=key,genomes=len(members),required=t,kmers=len(rep)))
        for s in rep: owners[s].add(key)
    for key,members in sorted(groups.items()):
        rep=core.sequence_dictionary(out/f'representative_kmers/{key}.fasta')
        sp=members[0]['species_taxid']
        extra=any(g['extra_species_background'] for g in members)
        allowed={s:h for s,h in rep.items() if len(owners[s])==1 and
                 (not extra or not (sp_owners.get(s,set())-{sp}))}
        rows.append(specificity_outputs(out,prep,members,allowed,key,'subgroup',members[0]['subgroup_label'],members[0]['taxid']))
    write_tsv(out/'representative_thresholds.tsv',['target_id','genomes','required','kmers'],threshold_rows)
    finish_bundle(out,rows)


def node_generate(args):
    prep,data=prepared(args.prepared)
    bynode=defaultdict(list)
    for g in data['genomes']:
        if args.node_rank not in g['lineage']:
            raise core.PipelineError(f"{g['accession']} has no {args.node_rank} ancestor")
        bynode[g['lineage'][args.node_rank]].append(g)
    selections={}
    if args.selection:
        for row in read_tsv(args.selection):
            node=int(row['node_taxid']); acc=row['accession_id'].strip()
            if node in selections: raise core.PipelineError(f'Duplicate manual node: {node}')
            if node not in bynode or acc not in {g['accession'] for g in bynode[node]}:
                raise core.PipelineError(f'Manual choice is not an eligible genome in node {node}: {acc}')
            selections[node]=acc
        missing=set(bynode)-set(selections)
        if missing and not args.allow_auto_unlisted:
            raise core.PipelineError(f'Manual selection must cover every node; missing: {sorted(missing)}')
    out=new_output(args); rows=[]; audit=[]
    for node,members in sorted(bynode.items()):
        g=next((g for g in members if g['accession']==selections.get(node)),None) or min(members,key=reference_key)
        seq=core.read_single_fasta(prep/g['normalized_fasta'])[1]
        length=args.window_length
        if len(seq)<length: raise core.PipelineError(f"Genome shorter than node window: {g['accession']}")
        key=f'{core.safe_component(args.node_rank)}_{node}'
        file=Path('probe_sets')/f'{key}.fasta'
        count=core.write_fasta(out/file,((window_header(g,s,length,'node_probe'),w)
                                      for s,w in circular_windows(seq,length,args.node_step)))
        rows.append(bundle_row(key,args.node_rank,data['names'][str(node)],node,file,count))
        audit.append(dict(node_taxid=node,accession_id=g['accession'],manual=node in selections,
                          non_atcg_fraction=g['non_atcg_fraction']))
    write_tsv(out/'selected_genomes.tsv',['node_taxid','accession_id','manual','non_atcg_fraction'],audit)
    finish_bundle(out,rows)


def input_targets(paths):
    result=[]
    for path in paths:
        path=Path(path).resolve()
        if path.is_dir():
            if not (path/'COMPLETE').is_file(): raise core.PipelineError(f'Incomplete bundle: {path}')
            meta=json.loads((path/'bundle.json').read_text())
            if meta['schema']!='emito.probes.v1': raise core.PipelineError(f'Unknown bundle: {path}')
            rows=read_tsv(path/'manifest.tsv')
            for row in rows:
                row=dict(row); row['path']=path/row['fasta']; result.append(row)
        else:
            if not path.is_file(): raise core.PipelineError(f'Input not found: {path}')
            result.append(dict(target_id=path.stem,rank='custom',name=path.stem,taxid='',path=path))
    for row in result:
        if not row['path'].is_file(): raise core.PipelineError(f"Missing target FASTA: {row['path']}")
    return result


def validate_probes(records):
    for h,s in records:
        if not core.DNA_RE.fullmatch(s):
            raise core.PipelineError(f'Probe is empty or contains non-ATCG bases: {h}')
        yield h,s


def collapse_records(records):
    grouped=defaultdict(list)
    for order,(header,seq) in enumerate(validate_probes(records)):
        acc=header.split('|',1)[0]
        fields=dict(re.findall(r'\|([A-Za-z_]+)=([^|\s]+)',header))
        try:
            start=int(fields.get('normalized_start',fields['start']))-1
            end=int(fields.get('normalized_end',fields['end']))-1
            n=int(fields['genome_length']) if 'genome_length' in fields else None
        except (KeyError,ValueError):
            raise core.PipelineError(f'Cannot parse collapse coordinates: {header}')
        if start<0 or end<start or end-start+1!=len(seq) or (n is not None and (start>=n or len(seq)>n)):
            raise core.PipelineError(f'Invalid/unrecognized coordinate interval: {header}')
        frame=fields.get('coordinate_frame','legacy')
        grouped[acc].append((start,end,order,header,seq,n,frame))
    kept=[]
    for acc,rows in sorted(grouped.items()):
        if len({(r[5],r[6]) for r in rows})!=1:
            raise core.PipelineError(f'Incompatible coordinate frames for {acc}; use the same prepare output')
        occupied=[]
        for start,end,order,h,s,n,frame in sorted(rows):
            parts=[(start,end)] if n is None or end<n else [(start,n-1),(0,end%n)]
            if any(a<=d and c<=b for a,b in parts for c,d in occupied): continue
            kept.append((h,s)); occupied.extend(parts)
    return kept


def access_or_collapse(args):
    targets=input_targets(args.inputs)
    if args.command=='access':
        if not 0<=args.gc_min<=args.gc_max<=100 or args.complexity_min>args.complexity_max:
            raise core.PipelineError('Invalid filter bounds')
        params=core.AccessParameters(args.gc_min,args.gc_max,args.complexity_min,args.complexity_max,args.dimer_k,args.dimer)
    out=new_output(args); rows=[]; counts=[]
    for i,row in enumerate(targets):
        key=f"{i+1:05d}_{core.safe_component(row['target_id'])}"
        file=Path('probe_sets')/f'{key}.{args.command}.fasta'
        if args.command=='access':
            stat=core.run_access_task(core.AccessTask(str(row['path']),str(out/file),str((out/file).with_suffix('.tsv')),params))
            count=stat['filtered']; incount=stat['input']
        else:
            records=list(core.iter_fasta(row['path']))
            count=core.write_fasta(out/file,collapse_records(records)); incount=len(records)
        rows.append(bundle_row(row['target_id'],row['rank'],row['name'],row['taxid'],file,count))
        counts.append(dict(input_fasta=row['path'],input_records=incount,output_records=count))
    write_tsv(out/'processing_summary.tsv',['input_fasta','input_records','output_records'],counts)
    finish_bundle(out,rows)


def merge(args):
    targets=input_targets(args.inputs)
    out=new_output(args); seen=set(); total=0; retained=0
    file=Path('probe_sets/merged.fasta'); (out/file).parent.mkdir()
    # Do not accumulate complete pools when deduplication is disabled.
    with (out/file).open('w') as f:
        for row in targets:
            for h,s in validate_probes(core.iter_fasta(row['path'])):
                total+=1
                if args.dedup and s in seen: continue
                if args.dedup: seen.add(s)
                f.write(f'>{h}\n{s}\n'); retained+=1
    write_tsv(out/'merge_inputs.tsv',['target_id','rank','name','taxid','path'],targets)
    json_write(out/'merge_summary.json',dict(input_records=total,output_records=retained,
        deduplicated=args.dedup,removed=total-retained))
    finish_bundle(out,[bundle_row('merged','merged','Merged probe pool','',file,retained)])


def positive(text):
    value=int(text)
    if value<1: raise argparse.ArgumentTypeError('Must be >=1')
    return value


def nonnegative(text):
    value=int(text)
    if value<0: raise argparse.ArgumentTypeError('Must be >=0')
    return value


def parser():
    p=argparse.ArgumentParser(description='eMito modular probe workflow')
    p.add_argument('--version',action='version',version=__version__)
    sub=p.add_subparsers(dest='command',required=True)
    prepare_p=sub.add_parser('prepare',help='QC, routing, normalization, optional alignment and circular windows')
    for opt in ['metadata','fasta-dir','names-dmp','nodes-dmp']:
        prepare_p.add_argument('--'+opt,type=Path,required=True)
    prepare_p.add_argument('--min-genome-length',type=positive,default=10000,
                          help='Minimum input genome length in bp (default: 10000); use a positive integer to adjust')
    prepare_p.add_argument('--align',action=argparse.BooleanOptionalAction,default=True,
                          help='Default: normalize strand/origin and align; --no-align skips ALL three')
    prepare_p.add_argument('--mafft',default='mafft')
    prepare_p.add_argument('--window-length',type=positive,default=52)
    prepare_p.add_argument('--probe-step',type=positive,default=5)
    prepare_p.add_argument('--kmer-step',type=positive,default=1)
    prepare_p.add_argument('--threads',type=positive,default=4)
    prepare_p.add_argument('--circular-anchor-length',type=positive,default=31)
    prepare_p.add_argument('--circular-min-anchor-length',type=positive,default=11)
    prepare_p.add_argument('--circular-min-anchor-hits',type=positive,default=3)
    taxa_p=sub.add_parser('taxa-generate',help='Species/genus plus subgroup specificity')
    taxa_p.add_argument('--prepared',type=Path,required=True)
    taxa_p.add_argument('--representative-fraction',type=float,default=.75)
    taxa_p.add_argument('--small-group-all-max',type=nonnegative,default=3)
    node_p=sub.add_parser('node-generate',help='Node-based circular tiling')
    node_p.add_argument('--prepared',type=Path,required=True)
    node_p.add_argument('--node-rank',default='family')
    node_p.add_argument('--selection',type=Path,help='TSV: node_taxid, accession_id; overrides automatic reference ranking')
    node_p.add_argument('--allow-auto-unlisted',action='store_true',help='Permit automatic selection for nodes absent from manual TSV')
    node_p.add_argument('--window-length',type=positive,default=52)
    node_p.add_argument('--node-step',type=positive,default=5)
    access_p=sub.add_parser('access',help='Assess each input target/pool independently')
    access_p.add_argument('--gc-min',type=float,default=35)
    access_p.add_argument('--gc-max',type=float,default=65)
    access_p.add_argument('--complexity-min',type=float,default=0)
    access_p.add_argument('--complexity-max',type=float,default=2)
    access_p.add_argument('--dimer-k',type=positive,default=11)
    access_p.add_argument('--dimer','--dimer-threshold',type=float,default=.15)
    collapse_p=sub.add_parser('collapse',help='Coordinate-greedy non-overlap selection, circular-aware')
    merge_p=sub.add_parser('merge',help='Merge target pools; exact-sequence deduplication by default')
    merge_p.add_argument('--dedup',action=argparse.BooleanOptionalAction,default=True)
    for child in [access_p,collapse_p,merge_p]:
        child.add_argument('--inputs',nargs='+',type=Path,required=True,help='Completed probe bundle directories or individual FASTAs')
    for child in [prepare_p,taxa_p,node_p,access_p,collapse_p,merge_p]:
        child.add_argument('--output',type=Path,required=True,help='NEW output directory; never silently overwritten')
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == 'legacy':
        from .legacy_cli import main as legacy_main
        return legacy_main(argv[1:])
    try:
        args=parser().parse_args(argv or ['--help'])
    except SystemExit as exc:
        return int(exc.code)
    try:
        for key,value in vars(args).items():
            if isinstance(value,float) and not math.isfinite(value):
                raise core.PipelineError(f'Non-finite value: {key}')
        dispatch={'prepare':prepare,'taxa-generate':taxa_generate,'node-generate':node_generate,
                  'access':access_or_collapse,'collapse':access_or_collapse,'merge':merge}
        dispatch[args.command](args)
        print(f'Completed {args.command}: {args.output.resolve()}')
        return 0
    except (core.PipelineError,OSError,ValueError,KeyError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        return 1


if __name__=='__main__':
    raise SystemExit(main())
