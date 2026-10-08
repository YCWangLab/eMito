import csv
import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from emito import cli


class ModularTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.root=Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def call(self,*args):
        self.assertEqual(cli.main(list(map(str,args))),0)

    def fixture(self):
        r=self.root
        taxa=[(1,1,'no rank','root'),(2,1,'family','Family'),(3,2,'genus','Genus'),
              (10,3,'species','Genus alpha'),(11,10,'subspecies','Genus alpha sub'),
              (20,3,'species','Genus beta'),(21,20,'subspecies','Genus beta sub'),
              (30,3,'species','Genus gamma'),(31,30,'subspecies','Genus gamma sub1'),
              (32,30,'subspecies','Genus gamma sub2'),(40,3,'species','Genus delta')]
        (r/'nodes.dmp').write_text(''.join(f'{t}\t|\t{p}\t|\t{rank}\t|\n' for t,p,rank,n in taxa))
        (r/'names.dmp').write_text(''.join(f'{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n' for t,p,rank,n in taxa))
        rng=random.Random(41)
        seq=''.join(rng.choice('ACGT') for _ in range(100))
        # Direct + lone subgroup, singleton subgroup fallback, multiple subgroup,
        # and explicit species-rank group exception.
        records=[('NC_000010.1','Genus alpha',10,'',seq),
                 ('AB000010.1','Genus alpha',10,'',seq[23:]+seq[:23]),
                 ('AB000011.1','Genus alpha',11,'alpha sub',seq),
                 ('AB000021.1','Genus beta',21,'beta sub','T'*100),
                 ('AB000031.1','Genus gamma',31,'gamma sub1','G'*100),
                 ('AB000032.1','Genus gamma',32,'gamma sub2','C'*100),
                 ('AB000040.1','Genus delta',40,'explicit population','A'*100)]
        (r/'fasta').mkdir()
        with (r/'metadata.tsv').open('w') as f:
            f.write('accession_id\tspecies_name\ttaxid\tsubgroup_label\n')
            for a,n,t,l,s in records:
                f.write(f'{a}\t{n}\t{t}\t{l}\n'); (r/'fasta'/f'{a}.fasta').write_text(f'>{a}\n{s}\n')
        return records

    def prep_args(self,output,align=False):
        return ['prepare','--metadata',self.root/'metadata.tsv','--fasta-dir',self.root/'fasta',
                '--names-dmp',self.root/'names.dmp','--nodes-dmp',self.root/'nodes.dmp',
                '--output',output,'--min-genome-length','50','--window-length','12',
                '--probe-step','5','--threads','1']+([] if align else ['--no-align'])

    def test_adjustable_minimum_genome_length(self):
        self.fixture()
        args = self.prep_args(self.root/'at_boundary')
        index = args.index('--min-genome-length') + 1
        args[index] = '100'
        self.call(*args)
        _, data = cli.prepared(self.root/'at_boundary')
        self.assertEqual(len(data['genomes']), 7)
        args[args.index('--output') + 1] = self.root/'too_short'
        args[index] = '101'
        self.assertEqual(cli.main(list(map(str, args))), 1)
        self.assertFalse((self.root/'too_short'/'COMPLETE').exists())
        for value in ('0', '-1'):
            args[index] = value
            self.assertEqual(cli.main(list(map(str, args))), 2)
        args[index] = '100'
        # Reusing preparation must not silently accept a new QC threshold.
        script = Path(__file__).resolve().parents[1]/'scripts/run_matched_comparison.py'
        result = subprocess.run([sys.executable, str(script), '--existing-standard',
            str(self.root), '--min-genome-length', '80', '--output', str(self.root/'compare')],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn('Cannot change the length filter', result.stderr)

    def records(self,path):
        return list(cli.core.iter_fasta(path))

    def test_defaults_threshold_and_reference(self):
        p=cli.parser()
        self.assertTrue(p.parse_args(['merge','--inputs','x','--output','y']).dedup)
        self.assertFalse(p.parse_args(['merge','--inputs','x','--output','y','--no-dedup']).dedup)
        for n,expected in [(1,1),(3,3),(4,3),(5,4),(6,5),(8,6)]:
            self.assertEqual(cli.core.representative_threshold(n,.75,3),expected)
        rows=[dict(accession='AB123.1',non_atcg_fraction=0),dict(accession='NC_123.1',non_atcg_fraction=.1),
              dict(accession='NC_124.1',non_atcg_fraction=.01)]
        self.assertEqual(min(rows,key=cli.reference_key)['accession'],'NC_124.1')
        self.assertEqual(min([dict(accession='AB123.1',non_atcg_fraction=.1),
                              dict(accession='AB124.1',non_atcg_fraction=0)],key=cli.reference_key)['accession'],'AB124.1')

    def test_no_align_raw_route_and_composition(self):
        original=self.fixture(); prep=self.root/'prep'
        self.call(*self.prep_args(prep))
        _,data=cli.prepared(prep)
        g={g['accession']:g for g in data['genomes']}
        for acc,n,t,label,seq in original:
            self.assertEqual(cli.core.read_single_fasta(prep/g[acc]['normalized_fasta'])[1],seq)
            self.assertEqual(g[acc]['rotation'],0)
            self.assertEqual(g[acc]['coordinate_frame'],'raw')
            self.assertEqual(g[acc]['probe_count'],20)  # includes end-to-start windows
        self.assertTrue(g['AB000011.1']['group_eligible'])
        self.assertTrue(g['AB000011.1']['extra_species_background'])
        self.assertTrue(g['AB000021.1']['taxa_eligible'])
        self.assertFalse(g['AB000021.1']['group_eligible'])
        self.assertFalse(g['AB000031.1']['taxa_eligible'])
        self.assertTrue(g['AB000031.1']['group_eligible'])
        self.assertFalse(g['AB000040.1']['taxa_eligible'])
        self.assertTrue(g['AB000040.1']['group_eligible'])
        taxa=self.root/'taxa'; node=self.root/'node'
        self.call('taxa-generate','--prepared',prep,'--output',taxa)
        self.call('node-generate','--prepared',prep,'--window-length',12,'--output',node)
        self.assertEqual(cli.read_tsv(node/'selected_genomes.tsv')[0]['accession_id'],'NC_000010.1')
        # Manual override can pick a non-NC record.
        (self.root/'choose.tsv').write_text('node_taxid\taccession_id\n2\tAB000021.1\n')
        self.call('node-generate','--prepared',prep,'--window-length',12,'--selection',self.root/'choose.tsv',
                  '--output',self.root/'manual')
        self.assertEqual(cli.read_tsv(self.root/'manual/selected_genomes.tsv')[0]['accession_id'],'AB000021.1')
        merged=self.root/'merge'; merged_raw=self.root/'merge_raw'
        self.call('merge','--inputs',taxa,node,'--output',merged)
        self.call('merge','--inputs',taxa,node,'--output',merged_raw,'--no-dedup')
        a=json.loads((merged/'merge_summary.json').read_text()); b=json.loads((merged_raw/'merge_summary.json').read_text())
        self.assertGreater(b['output_records'],a['output_records'])
        self.call('access','--inputs',merged,'--output',self.root/'access','--gc-min',0,'--gc-max',100,
                  '--complexity-max',100,'--dimer',0)
        self.call('collapse','--inputs',self.root/'access','--output',self.root/'collapsed')
        self.call('collapse','--inputs',taxa,'--output',self.root/'collapse_first')
        self.call('access','--inputs',self.root/'collapse_first','--output',self.root/'access_second')
        self.call('merge','--inputs',self.root/'access_second',node,'--output',self.root/'final')
        self.assertTrue((self.root/'final/COMPLETE').is_file())

    def test_special_subgroup_other_species_background(self):
        # Build minimal prepared k-mer/probe sets to isolate membership semantics.
        p=self.root/'p'; p.mkdir(); (p/'COMPLETE').write_text('prepare')
        genomes=[]
        specs=[('AB001.1',10,10,'',True,False,False,['AAAA','CCCC','GGGG']),
               ('AB002.1',20,20,'',True,False,False,['AAAA']),
               ('AB003.1',10,11,'only subgroup',False,True,True,['AAAA','CCCC','GGGG']),
               ('AB004.1',30,31,'other subgroup',False,True,False,['GGGG'])]
        for acc,sp,tid,label,ta,gr,extra,seqs in specs:
            path=f'{acc}.fa'; cli.core.write_fasta(p/path,[(acc,s) for s in seqs])
            genomes.append(dict(accession=acc,species_taxid=sp,taxid=tid,genus_taxid=3,species_name=f'Species {sp}',
                subgroup_label=label,taxa_eligible=ta,group_eligible=gr,extra_species_background=extra,
                group_target=f'group_{tid}',probe_file=path,kmer_file=path))
        cli.json_write(p/'prepared.json',dict(schema='emito.prepared.v1',genomes=genomes))
        out=self.root/'out'; self.call('taxa-generate','--prepared',p,'--output',out)
        # AAAA excluded by non-focal species, GGGG excluded by other subgroup,
        # CCCC permitted even though also present in direct focal species.
        self.assertEqual([s for _,s in self.records(out/'probe_sets/group_11.fasta')],['CCCC'])

    def test_circular_window_and_collapse(self):
        self.assertEqual(list(cli.circular_windows('ACGTAC',4,5)),[(0,'ACGT'),(5,'CACG')])
        self.assertFalse(any('N' in s for _,s in cli.circular_windows('ACGNAC',4,1)))
        g=dict(accession='AB001.1',length=100,orientation='Forward',rotation=0,coordinate_frame='raw')
        # Starts 0, 49, 52, 90: 0..51 retained; last wraps and overlaps first.
        records=[(cli.window_header(g,s,52,'probe'),'A'*52) for s in [0,49,52,90]]
        self.assertEqual(len(cli.collapse_records(records)),1)
        records=[(cli.window_header(g,s,20,'probe'),'A'*20) for s in [0,20,40,60,80]]
        self.assertEqual(len(cli.collapse_records(records)),5)
        other=dict(g,accession='AB002.1')
        self.assertEqual(len(cli.collapse_records([records[0],(cli.window_header(other,0,20,'probe'),'A'*20)])),2)
        self.assertEqual(len(cli.collapse_records([records[-1],records[-1]])),1)

    def test_alignment_normalization_and_wrap(self):
        # Mock only external MAFFT: identical normalized sequences need no gaps.
        root=self.root
        seq=''.join(random.Random(5+i).choice('ACGT') for i in range(100))
        items=[]
        for acc,s in [('NC_001.1',seq),('AB002.1',seq[17:]+seq[:17])]:
            p=root/f'{acc}.fa'; p.write_text(f'>{acc}\n{s}\n')
            items.append(dict(accession=acc,organized_fasta=str(p),alignment_group='test',length=100,
                              non_atcg_fraction=0))
        def fake(cmd,stdout,stderr):
            stdout.write(Path(cmd[-1]).read_text())
            return subprocess.CompletedProcess(cmd,0)
        with patch.object(cli.subprocess,'run',side_effect=fake):
            members=cli.prepare_group((items,str(root),True,'mafft',15,7,3,12,5,1))
        normalized=[cli.core.read_single_fasta(root/g['normalized_fasta'])[1] for g in members]
        self.assertEqual(normalized,[seq,seq])
        self.assertTrue(any('|wrap=1|' in h for h,s in self.records(root/members[0]['probe_file'])))
        self.assertEqual([s for h,s in self.records(root/members[0]['probe_file'])],
                         [s for h,s in self.records(root/members[1]['probe_file'])])

    def test_ambiguous_genome_is_preserved(self):
        seq = 'ACGT' * 30
        ambiguous = seq[:40] + 'N' + seq[41:]
        for align in (False, True):
            root = self.root / str(align)
            root.mkdir()
            items = []
            for acc, s in [('NC_001.1',seq), ('AB002.1',ambiguous)]:
                p = root / (acc + '.fa')
                p.write_text(f'>{acc}\n{s}\n')
                items.append(dict(accession=acc,organized_fasta=str(p),alignment_group='test',
                                  length=len(s),non_atcg_fraction=cli.fraction_bad(s)))
            def fake(cmd,stdout,stderr):
                stdout.write(Path(cmd[-1]).read_text())
                return subprocess.CompletedProcess(cmd,0)
            with patch.object(cli.subprocess,'run',side_effect=fake), \
                 patch.object(cli.core,'circular_normalize_to_reference',side_effect=lambda r,s,*a:
                    type('Result',(),dict(sequence=s,orientation='Forward',rotation=0,
                                         anchor_length=11,anchor_hits=10))()):
                members = cli.prepare_group((items,str(root),align,'mafft',31,11,3,12,5,1))
            self.assertEqual(len(self.records(root/'alignments/test/normalized_input.fasta')),2)
            self.assertEqual(cli.core.read_single_fasta(root/'normalized/AB002.1.fasta')[1],ambiguous)
            probes = self.records(root/members[1]['probe_file'])
            self.assertTrue(probes)
            self.assertTrue(all(set(s) <= set('ACGT') for h,s in probes))

    def test_refseq_audit_and_matched_runner(self):
        records = self.fixture()
        scripts = Path(__file__).resolve().parents[1]/'scripts'
        fna = self.root/'all.fna'
        fna.write_text(''.join(f'>{a}\n{s}\n' for a,n,t,l,s in records))
        mapping = self.root/'mapping.tsv'
        mapping.write_text('accession\toriginal_taxid\n'+''.join(f'{a}\t{t}\n' for a,n,t,l,s in records))
        audit = self.root/'audit'
        subprocess.run([sys.executable,str(scripts/'audit_prepare_refseq.py'),
            '--fasta',str(fna),'--mapping',str(mapping),'--nodes-dmp',str(self.root/'nodes.dmp'),
            '--names-dmp',str(self.root/'names.dmp'),'--output',str(audit),'--prepare',
            '--min-genome-length','50'],check=True,capture_output=True)
        stats = {r['metric']:r['value'] for r in cli.read_tsv(audit/'summary.tsv')}
        self.assertEqual(stats['species_with_multiple_records'],'2')
        self.assertEqual(stats['eligible_records'],'7')
        std=self.root/'std'; std.mkdir()
        self.call(*self.prep_args(std/'01_prepare'))
        self.call('taxa-generate','--prepared',std/'01_prepare','--output',std/'02_taxa_generate')
        self.call('access','--inputs',std/'02_taxa_generate','--output',std/'03_taxa_access')
        result=self.root/'matched'
        subprocess.run([sys.executable,str(scripts/'run_matched_comparison.py'),
                        '--existing-standard',str(std),'--output',str(result)],check=True,capture_output=True)
        self.assertTrue((result/'COMPLETE').is_file())
        self.assertTrue(all(r['rank']=='species' for r in cli.read_tsv(result/'03_taxa_access/manifest.tsv')))
        self.assertEqual(len(cli.read_tsv(result/'03_taxa_access/manifest.tsv')),2)
        single=self.root/'single.tsv'
        single.write_text('accession_id\tspecies_name\ttaxid\tsubgroup_label\n'
                          'NC_000010.1\tGenus alpha\t10\t\n'
                          'AB000021.1\tGenus beta\t21\t\n')
        singleout=self.root/'single_comparison'
        args=[sys.executable,str(scripts/'run_matched_comparison.py'),
              '--metadata',str(single),'--fasta-dir',str(self.root/'fasta'),
              '--nodes-dmp',str(self.root/'nodes.dmp'),'--names-dmp',str(self.root/'names.dmp'),
              '--no-align','--require-single-reference','--min-genome-length','1','--threads','1',
              '--mafft','does-not-exist','--output',str(singleout)]
        subprocess.run(args,check=True,capture_output=True)
        self.assertTrue((singleout/'COMPLETE').is_file())
        data=json.loads((singleout/'01_prepare/prepared.json').read_text())
        self.assertFalse(data['aligned'])
        self.assertEqual(len(data['genomes']),2)
        with single.open('a') as f:
            f.write('AB000010.1\tGenus alpha\t10\t\n')
        args[-1]=str(self.root/'invalid_single')
        failed=subprocess.run(args,capture_output=True,text=True)
        self.assertNotEqual(failed.returncode,0)
        self.assertIn('Expected one reference per species',failed.stderr)

    def test_no_align_never_normalizes(self):
        items=[]
        for acc,seq in [('NC_001.1','A'*80),('AB002.1','C'*80)]:
            p=self.root/f'{acc}.fa'; p.write_text(f'>{acc}\n{seq}\n')
            items.append(dict(accession=acc,organized_fasta=str(p),alignment_group='test',length=80,non_atcg_fraction=0))
        with patch.object(cli.core,'circular_normalize_to_reference',side_effect=AssertionError('called normalization')), \
             patch.object(cli.subprocess,'run',side_effect=AssertionError('called MAFFT')):
            cli.prepare_group((items,str(self.root),False,None,31,11,3,52,5,1))


if __name__=='__main__':
    unittest.main()
