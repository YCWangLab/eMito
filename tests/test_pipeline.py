#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from emito import access_reporting
from emito import pipeline
from emito import cli
from emito import reporting as mode_merge_statistics
from emito import stage_reporting as probe_statistics
from emito import taxonomy_reporting as input_taxonomy_statistics


SCRIPT = Path(pipeline.__file__)


class EMitoPipelineTest(unittest.TestCase):
    def write(self, path: Path, text: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def add_taxon(
        self,
        nodes: list[str],
        names: list[str],
        taxid: int,
        parent: int,
        rank: str,
        name: str,
    ) -> None:
        nodes.append(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\n")
        names.append(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\t\n")

    def test_representative_thresholds(self) -> None:
        expected = {1: 1, 2: 2, 3: 3, 4: 3, 5: 4, 6: 5, 8: 6}
        observed = {
            count: pipeline.representative_threshold(count, 0.75, 3)
            for count in expected
        }
        self.assertEqual(observed, expected)

    def test_public_cli_version_and_unknown_command(self) -> None:
        self.assertEqual(cli.main(["--version"]), 0)
        self.assertEqual(cli.main(["not-a-command"]), 2)

    def test_circular_normalization_handles_rotation_and_reverse_complement(self) -> None:
        reference = (
            "ACGTTGCATGTCAGTACCGATGACCTAGGTCATCGTACGATTCGAGCCTATGCACTGATC"
            "GGTACCTTAGCGATACGTCAGGATCCGTTACGATGCTAGTCGATCGTACCTGACTGACGA"
        )
        rotation = 37
        rotated = reference[rotation:] + reference[:rotation]
        result = pipeline.circular_normalize_to_reference(
            reference, rotated, 13, 9, 3
        )
        self.assertEqual(result.sequence, reference)
        self.assertEqual(result.orientation, "Forward")
        self.assertGreaterEqual(result.modal_anchor_support, 3)

        reverse_rotated = pipeline.reverse_complement(rotated)
        result = pipeline.circular_normalize_to_reference(
            reference, reverse_rotated, 13, 9, 3
        )
        self.assertEqual(result.sequence, reference)
        self.assertEqual(result.orientation, "ReverseComplement")
        self.assertGreaterEqual(result.modal_anchor_support, 3)

    def test_circular_normalization_adaptively_falls_back_to_shorter_anchors(self) -> None:
        reference = (
            "ACGTTGCATGTCAGTACCGATGACCTAGGTCATCGTACGATTCGAGCCTATGCACTGATC"
            "GGTACCTTAGCGATACGTCAGGATCCGTTACGATGCTAGTCGATCGTACCTGACTGACGA"
        )
        # Mutations every 19 bases disrupt every sampled 31-mer while leaving
        # enough shorter exact anchors for reliable orientation/origin voting.
        mutated = list(reference)
        for index in range(9, len(mutated), 19):
            mutated[index] = {"A": "C", "C": "G", "G": "T", "T": "A"}[mutated[index]]
        mutated_sequence = "".join(mutated)
        rotation = 29
        rotated = mutated_sequence[rotation:] + mutated_sequence[:rotation]
        result = pipeline.circular_normalize_to_reference(
            reference, rotated, 31, 11, 3
        )
        self.assertLess(result.anchor_length, 31)
        self.assertEqual(result.sequence, mutated_sequence)
        self.assertEqual(result.orientation, "Forward")
        self.assertEqual(result.attempted_hits[0][0], 31)
        self.assertGreaterEqual(result.anchor_hits, 3)

    def test_alignment_coordinate_windows_do_not_follow_indel_shifted_phase(self) -> None:
        # Target has one inserted G after reference position 3. Shared reference
        # starts at offsets 0,2,4 map to target offsets 0,2,5 (not 0,2,4).
        reference_alignment = "ACG-TACGTACG"
        target_alignment = "ACGGTACGTACG"
        reference_columns = pipeline.aligned_reference_columns(reference_alignment)
        ungapped, column_to_offset, base_columns = pipeline.alignment_coordinate_maps(
            target_alignment
        )
        target_offsets = []
        for reference_offset in range(0, 7, 2):
            window = pipeline.extract_aligned_window(
                target_alignment,
                ungapped,
                column_to_offset,
                base_columns,
                reference_columns[reference_offset],
                3,
            )
            self.assertIsNotNone(window)
            target_offsets.append(window[1])
        self.assertEqual(target_offsets, [0, 2, 5, 7])

    def test_window_task_uses_mafft_reference_columns_across_an_indel(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            reference = (
                "ACGTTGCATGTCAGTACCGATGACCTAGGTCATCGTACGATTCGAGCCTATGCACTGATC"
                "GGTACCTTAGCGATACGTCAGGATCCGTTACGATGCTAGTCGATCGTACCTGACTGACGA"
            )
            insertion_offset = 45
            target = reference[:insertion_offset] + "A" + reference[insertion_offset:]
            reference_path = root / "NC_000001.1.fasta"
            target_path = root / "AB000002.1.fasta"
            self.write(reference_path, f">NC_000001.1\n{reference}\n")
            self.write(target_path, f">AB000002.1\n{target}\n")

            fake_mafft = root / "fake_mafft.py"
            self.write(
                fake_mafft,
                f"""#!{sys.executable}
import pathlib
import sys

records = {{}}
header = None
parts = []
for line in pathlib.Path(sys.argv[-1]).read_text(encoding="ascii").splitlines():
    if line.startswith(">"):
        if header is not None:
            records[header] = "".join(parts)
        header = line[1:]
        parts = []
    else:
        parts.append(line.strip())
records[header] = "".join(parts)
reference = records["NC_000001.1"]
target = records["AB000002.1"]
offset = {insertion_offset}
print(">NC_000001.1")
print(reference[:offset] + "-" + reference[offset:])
print(">AB000002.1")
print(target)
""",
            )
            fake_mafft.chmod(0o755)

            def member(accession: str, source: Path) -> pipeline.WindowMember:
                output = root / "windows" / accession
                return pipeline.WindowMember(
                    accession=accession,
                    input_fasta=str(source),
                    probe_output=str(output.with_suffix(".probe.fasta")),
                    kmer_forward_output=str(output.with_suffix(".kmer.forward.fasta")),
                    kmer_reverse_output=str(
                        output.with_suffix(".kmer.reverse_complement.fasta")
                    ),
                    kmer_merged_output=str(output.with_suffix(".kmer.merged.fasta")),
                )

            reference_member = member("NC_000001.1", reference_path)
            target_member = member("AB000002.1", target_path)
            task = pipeline.WindowTask(
                group_id="species_test",
                members=(reference_member, target_member),
                reference_accession="NC_000001.1",
                normalized_input=str(root / "alignment" / "normalized.fasta"),
                alignment_output=str(root / "alignment" / "aligned.fasta"),
                normalization_manifest=str(root / "alignment" / "normalization.tsv"),
                mafft_executable=str(fake_mafft),
                circular_anchor_length=13,
                circular_min_anchor_length=9,
                circular_min_anchor_hits=3,
                window_length=10,
                probe_step=5,
                kmer_step=1,
            )
            result = pipeline.run_window_task(task)
            self.assertEqual(result["genomes"], 2)

            reference_records = list(
                pipeline.iter_fasta(Path(reference_member.probe_output))
            )
            target_records = list(pipeline.iter_fasta(Path(target_member.probe_output))
            )
            self.assertEqual(len(reference_records), len(target_records))
            post_indel_index = next(
                index
                for index, (header, _) in enumerate(reference_records)
                if "|normalized_start=51|" in header
            )
            self.assertIn(
                "|normalized_start=52|", target_records[post_indel_index][0]
            )
            self.assertEqual(
                reference_records[post_indel_index][1],
                target_records[post_indel_index][1],
            )

    def test_collapse_matches_requested_greedy_example(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "source.fasta"
            output = root / "collapsed.fasta"
            self.write(
                source,
                """>ACC1|start=1|end=4|strand=Forward
AAAA
>ACC1|start=2|end=5|strand=Forward
AAAT
>ACC1|start=4|end=7|strand=Forward
AATT
>ACC1|start=5|end=8|strand=Forward
ATTT
>ACC2|start=1|end=4|strand=Forward
CCCC
""",
            )
            result = pipeline.run_collapse_task(
                pipeline.CollapseTask(str(source), str(output))
            )
            records = list(pipeline.iter_fasta(output))
            self.assertEqual(result, {"input": 5, "retained": 3})
            self.assertEqual([sequence for _, sequence in records], ["AAAA", "ATTT", "CCCC"])

    def test_final_merge_dedup_can_be_disabled(self) -> None:
        self.assertTrue(pipeline.parse_args([]).final_dedup)
        self.assertFalse(pipeline.parse_args(["--no-final-dedup"]).final_dedup)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            first = root / "first.fasta"
            second = root / "second.fasta"
            self.write(first, ">A|start=1|end=4\nAAAA\n>B|start=1|end=4\nCCCC\n")
            self.write(second, ">C|start=1|end=4\nAAAA\n")
            first_target = pipeline.Target("species_1", "species", "Species one", 1)
            second_target = pipeline.Target("species_2", "species", "Species two", 2)
            mode_outputs = {
                "taxa": [(first_target, first), (second_target, second)]
            }

            dedup_path, dedup_count = pipeline.merge_final_outputs(
                mode_outputs, root / "dedup", deduplicate=True
            )
            merged_path, merged_count = pipeline.merge_final_outputs(
                mode_outputs, root / "not_dedup", deduplicate=False
            )
            self.assertEqual(dedup_count, 2)
            self.assertEqual(merged_count, 3)
            self.assertEqual(len(list(pipeline.iter_fasta(dedup_path))), 2)
            self.assertEqual(len(list(pipeline.iter_fasta(merged_path))), 3)
            self.assertEqual(
                [sequence for _, sequence in pipeline.iter_fasta(merged_path)],
                ["AAAA", "CCCC", "AAAA"],
            )

    def test_eprobe_access_defaults_and_dimer_self_subtraction(self) -> None:
        args = pipeline.parse_args([])
        self.assertEqual((args.gc_min, args.gc_max), (35.0, 65.0))
        self.assertEqual((args.complexity_min, args.complexity_max), (0.0, 2.0))
        self.assertEqual(args.dimer_k, 11)
        self.assertEqual(args.dimer_threshold, 0.15)

        self.assertEqual(
            pipeline.dimer_scores([("probe_a", "AAAACCCC")], 4),
            [0.0],
        )
        paired = pipeline.dimer_scores(
            [("probe_a", "AAAACCCC"), ("probe_b", "GGGGTTTT")], 4
        )
        self.assertTrue(all(score > 0 for score in paired))
        self.assertEqual(pipeline.eprobe_dimer_cutoff([1.0, 2.0, 3.0, 4.0], 0.5), 3.0)
        self.assertEqual(pipeline.eprobe_dimer_cutoff([1.0, 2.0], 5.0), 5.0)

    def test_complete_pipeline(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            metadata = root / "metadata.tsv"
            fasta_dir = root / "fasta"
            nodes_path = root / "nodes.dmp"
            names_path = root / "names.dmp"
            output_root = root / "output"
            fake_mafft = root / "fake_mafft.py"
            self.write(
                fake_mafft,
                f"""#!{sys.executable}
import pathlib
import sys

# Test genomes in this end-to-end fixture have equal lengths, so passing the
# normalized input through is a valid alignment for exercising pipeline logic.
sys.stdout.write(pathlib.Path(sys.argv[-1]).read_text(encoding="ascii"))
""",
            )
            fake_mafft.chmod(0o755)

            nodes: list[str] = []
            names: list[str] = []
            self.add_taxon(nodes, names, 1, 1, "no rank", "root")
            self.add_taxon(nodes, names, 10, 1, "family", "Familia one")
            self.add_taxon(nodes, names, 20, 10, "genus", "GenusA")
            self.add_taxon(nodes, names, 21, 20, "species", "GenusA alpha")
            self.add_taxon(nodes, names, 22, 21, "subspecies", "GenusA alpha subA")
            self.add_taxon(nodes, names, 23, 21, "subspecies", "GenusA alpha subB")
            self.add_taxon(nodes, names, 24, 20, "species", "GenusA delta")
            self.add_taxon(nodes, names, 30, 10, "genus", "GenusB")
            self.add_taxon(nodes, names, 31, 30, "species", "GenusB beta")
            self.add_taxon(nodes, names, 40, 1, "family", "Familia two")
            self.add_taxon(nodes, names, 41, 40, "genus", "GenusC")
            self.add_taxon(nodes, names, 42, 41, "species", "GenusC gamma")
            self.add_taxon(nodes, names, 50, 40, "genus", "GenusD")
            self.add_taxon(nodes, names, 51, 50, "species", "GenusD epsilon")
            self.add_taxon(nodes, names, 52, 51, "subspecies", "GenusD epsilon subOnly")
            self.add_taxon(nodes, names, 60, 40, "genus", "GenusE")
            self.add_taxon(nodes, names, 61, 60, "species", "GenusE zeta")
            self.write(nodes_path, "".join(nodes))
            self.write(names_path, "".join(names))

            rows = [
                ("NC_000001.1", "GenusA alpha", 21, "", "AAAACCCCGGGG"),
                ("AB000002.1", "GenusA alpha", 21, "", "AAAACCCCTTTT"),
                ("AB000003.1", "GenusB beta", 31, "", "GATCGATCGATC"),
                ("AB000004.1", "GenusC gamma", 42, "", "TGCATGCATGCA"),
                # This second GenusA species lets the test distinguish genus-
                # specific k-mers from species-specific k-mers.
                ("AB000009.1", "GenusA delta", 24, "", "AAAACCCCAAAA"),
                ("AB000005.1", "GenusA alpha", 22, "GenusA alpha subA", "TTTTAAAACCCC"),
                ("AB000006.1", "GenusA alpha", 22, "GenusA alpha subA", "TTTTAAAACCCA"),
                ("AB000007.1", "GenusA alpha", 23, "GenusA alpha subB", "CGCGGGGGTTTT"),
                # N-containing windows must be dropped from all generated window files.
                ("AB000008.1", "GenusA alpha", 23, "GenusA alpha subB", "CGCGNNNNATAT"),
                # With no direct species genome and only one true taxonomic
                # subgroup, these records must be routed to taxa as a species
                # fallback and omitted from group mode.
                ("AB000011.1", "GenusD epsilon", 52, "GenusD epsilon subOnly", "ACGTACGTAAAA"),
                ("AB000012.1", "GenusD epsilon", 52, "GenusD epsilon subOnly", "ACGTACGTAAAT"),
                # A labelled record whose actual TaxID is already the species
                # TaxID is an explicit group/population, not a taxonomic
                # subgroup fallback. This models the Heidelbergensis exception.
                ("AB000013.1", "GenusE zeta", 61, "Special population", "GGGGAAAATTTT"),
                # Input-length QC must remove this partial record before it is
                # organized, aligned or counted as an individual.
                ("AB000010.1", "GenusA alpha", 21, "", "AA"),
            ]
            metadata.parent.mkdir(parents=True, exist_ok=True)
            with metadata.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(("accession_id", "species_name", "taxid", "subgroup_label"))
                for accession, species, taxid, subgroup, sequence in rows:
                    writer.writerow((accession, species, taxid, subgroup))
                    self.write(fasta_dir / f"{accession}.fasta", f">source {accession}\n{sequence}\n")

            command = [
                sys.executable,
                str(SCRIPT),
                "--metadata",
                str(metadata),
                "--fasta-dir",
                str(fasta_dir),
                "--nodes-dmp",
                str(nodes_path),
                "--names-dmp",
                str(names_path),
                "--output-root",
                str(output_root),
                "--window-length",
                "4",
                "--probe-step",
                "2",
                "--kmer-step",
                "1",
                "--min-genome-length",
                "4",
                "--mafft",
                str(fake_mafft),
                "--circular-anchor-length",
                "4",
                "--circular-min-anchor-length",
                "4",
                "--circular-min-anchor-hits",
                "1",
                "--node-step",
                "2",
                "--threads",
                "1",
                "--gc-max",
                "100",
                "--complexity-max",
                "100",
                "--dimer-k",
                "2",
                "--dimer-threshold",
                "0",
                "--group-collapse",
            ]
            completed = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)

            organized = (
                output_root
                / "01_organized_genomes"
                / "GenusA"
                / "GenusA_alpha"
                / "GenusA_alpha_subA"
                / "AB000005.1.fasta"
            )
            self.assertTrue(organized.is_file())
            self.assertFalse(
                any(
                    (output_root / "01_organized_genomes").rglob(
                        "AB000010.1.fasta"
                    )
                )
            )
            with (output_root / "00_input_qc" / "genome_length_qc.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                length_qc = {
                    row["accession_id"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(length_qc["AB000010.1"]["sequence_length"], "2")
            self.assertEqual(length_qc["AB000010.1"]["status"], "excluded")
            excluded_qc = (
                output_root / "00_input_qc" / "excluded_short_genomes.tsv"
            ).read_text(encoding="utf-8")
            self.assertIn("AB000010.1", excluded_qc)

            reverse_file = next(
                (output_root / "02_windows").rglob("NC_000001.1.kmer.reverse_complement.fasta")
            )
            reverse_records = list(pipeline.iter_fasta(reverse_file))
            self.assertTrue(reverse_records)
            self.assertTrue(all("[Reverse Complement]" in header for header, _ in reverse_records))
            self.assertTrue(all("start=" in header and "end=" in header for header, _ in reverse_records))

            invalid_source_outputs = list(
                (output_root / "02_windows").rglob("AB000008.1.*.fasta")
            )
            self.assertEqual(len(invalid_source_outputs), 4)
            for path in invalid_source_outputs:
                self.assertTrue(
                    all("N" not in sequence for _, sequence in pipeline.iter_fasta(path)), path
                )

            selection = output_root / "05_eMito_node_generate" / "selected_genomes.tsv"
            with selection.open(encoding="utf-8", newline="") as handle:
                selected = list(csv.DictReader(handle, delimiter="\t"))
            family_one = next(row for row in selected if row["node_taxid"] == "10")
            self.assertEqual(family_one["selected_accession"], "NC_000001.1")

            taxa_access_files = list(
                (output_root / "06_optional_modes" / "taxa" / "access").glob("*.fasta")
            )
            group_collapse_files = list(
                (output_root / "06_optional_modes" / "group" / "collapse").glob("*.fasta")
            )
            self.assertTrue(taxa_access_files)
            self.assertFalse(
                any("__genus_taxid_" in path.name for path in taxa_access_files)
            )
            self.assertTrue(
                all("__species_taxid_" in path.name for path in taxa_access_files)
            )
            self.assertTrue(group_collapse_files)
            self.assertFalse((output_root / "06_optional_modes" / "node").exists())

            taxa_intersections = list(
                (output_root / "03_eMito_taxa_generate" / "intersections").rglob(
                    "*.intersection.probe.fasta"
                )
            )
            group_intersections = list(
                (output_root / "04_eMito_group_generate" / "intersections").rglob(
                    "*.intersection.probe.fasta"
                )
            )
            self.assertTrue(taxa_intersections)
            self.assertTrue(group_intersections)

            taxa_root = output_root / "03_eMito_taxa_generate"
            with (taxa_root / "probe_sets" / "manifest.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                taxa_manifest = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(taxa_manifest), 5)
            self.assertTrue(all(row["rank"] == "species" for row in taxa_manifest))
            self.assertTrue(
                all(row["target_id"].startswith("species_taxid_") for row in taxa_manifest)
            )
            self.assertFalse((taxa_root / "probe_sets" / "genus").exists())
            self.assertFalse((taxa_root / "intersections" / "genus").exists())

            genus_kmer_files = list(
                (taxa_root / "specific_kmers" / "genus").glob(
                    "*.genus_specific_kmer.fasta"
                )
            )
            combined_kmer_files = list(
                (taxa_root / "specific_kmers" / "combined_by_species").glob(
                    "*.taxa_specific_kmer.fasta"
                )
            )
            self.assertEqual(len(genus_kmer_files), 4)
            self.assertEqual(len(combined_kmer_files), 5)

            # A sequence shared by the two GenusA species is genus-specific but
            # not species-specific, so it must enter alpha's combined set.
            alpha_species_file = next(
                (taxa_root / "specific_kmers" / "species").glob(
                    "GenusA_alpha__*.species_specific_kmer.fasta"
                )
            )
            alpha_combined_file = next(
                (taxa_root / "specific_kmers" / "combined_by_species").glob(
                    "GenusA_alpha__*.taxa_specific_kmer.fasta"
                )
            )
            alpha_species_sequences = {
                sequence for _, sequence in pipeline.iter_fasta(alpha_species_file)
            }
            alpha_combined_sequences = {
                sequence for _, sequence in pipeline.iter_fasta(alpha_combined_file)
            }
            self.assertGreater(
                len(alpha_combined_sequences - alpha_species_sequences), 0
            )
            for combined_path in combined_kmer_files:
                species_component = combined_path.name.split(".taxa_specific_kmer", 1)[0]
                representative_path = (
                    taxa_root
                    / "representative_kmers"
                    / f"{species_component}.representative_kmer.fasta"
                )
                combined_sequences = {
                    sequence for _, sequence in pipeline.iter_fasta(combined_path)
                }
                representative_sequences = {
                    sequence for _, sequence in pipeline.iter_fasta(representative_path)
                }
                self.assertTrue(combined_sequences <= representative_sequences)

            # Taxa uses direct species genomes plus a lone true-taxonomic-
            # subgroup fallback. Multi-subgroup and explicit population inputs
            # remain excluded from taxa.
            comparative_subgroup_accessions = {
                "AB000005.1",
                "AB000006.1",
                "AB000007.1",
                "AB000008.1",
            }
            fallback_accessions = {"AB000011.1", "AB000012.1"}
            explicit_group_accessions = {"AB000013.1"}
            taxa_representatives = list(
                (output_root / "03_eMito_taxa_generate" / "representative_kmers").glob(
                    "*.representative_kmer.fasta"
                )
            )
            self.assertTrue(taxa_representatives)
            taxa_headers = {
                header.split("|", 1)[0]
                for path in taxa_representatives
                for header, _ in pipeline.iter_fasta(path)
            }
            self.assertTrue(taxa_headers.isdisjoint(comparative_subgroup_accessions))
            self.assertTrue(taxa_headers.isdisjoint(explicit_group_accessions))
            self.assertTrue(taxa_headers & fallback_accessions)

            group_representatives = list(
                (output_root / "04_eMito_group_generate" / "representative_kmers").rglob(
                    "*.representative_kmer.fasta"
                )
            )
            group_headers = {
                header.split("|", 1)[0]
                for path in group_representatives
                for header, _ in pipeline.iter_fasta(path)
            }
            self.assertTrue(group_headers & comparative_subgroup_accessions)
            self.assertTrue(group_headers & explicit_group_accessions)
            self.assertTrue(group_headers.isdisjoint(fallback_accessions))

            with (output_root / "00_input_qc" / "mode_routing.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                routing = {
                    row["species_name"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(
                routing["GenusD epsilon"]["taxa_input_decision"],
                "single_taxonomic_subgroup_fallback",
            )
            self.assertEqual(
                routing["GenusD epsilon"]["group_input_decision"],
                "single_taxonomic_subgroup_not_comparable",
            )
            self.assertEqual(
                routing["GenusE zeta"]["taxa_input_decision"],
                "excluded_explicit_group_without_direct_genome",
            )
            self.assertEqual(
                routing["GenusE zeta"]["group_input_decision"],
                "explicit_user_defined_group",
            )

            final_fasta = output_root / "07_final_probe_set" / "final_probe_set.fasta"
            final_records = list(pipeline.iter_fasta(final_fasta))
            self.assertTrue(final_records)
            self.assertEqual(
                len(final_records), len({sequence for _, sequence in final_records})
            )
            manifest_text = (
                output_root / "07_final_probe_set" / "input_probe_files.tsv"
            ).read_text(encoding="utf-8")
            self.assertIn(".assessed.fasta", manifest_text)
            self.assertIn(".collapsed.fasta", manifest_text)
            self.assertIn(".node_tiling.probe.fasta", manifest_text)

            report_dir = root / "statistics"
            self.assertEqual(
                probe_statistics.main(
                    [
                        "--output-root",
                        str(output_root),
                        "--report-dir",
                        str(report_dir),
                    ]
                ),
                0,
            )
            with (report_dir / "probe_stage_summary.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                summary = {
                    row["stage"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(
                summary["11_terminal_outputs_all_modes"][
                    "unique_after_merging_all_stage_files"
                ],
                summary["12_final_probe_set"][
                    "unique_after_merging_all_stage_files"
                ],
            )
            self.assertTrue((report_dir / "probe_file_counts.tsv").is_file())
            self.assertTrue((report_dir / "kmer_stage_summary.tsv").is_file())
            self.assertTrue((report_dir / "kmer_file_counts.tsv").is_file())
            self.assertTrue(
                (report_dir / "probe_terminal_mode_contributions.tsv").is_file()
            )

            mode_report = report_dir / "mode_processing_summary.tsv"
            final_merge_report = report_dir / "final_merge_summary.tsv"
            self.assertEqual(
                mode_merge_statistics.main(
                    [
                        "--output-root",
                        str(output_root),
                        "--report",
                        str(mode_report),
                        "--final-report",
                        str(final_merge_report),
                    ]
                ),
                0,
            )
            with mode_report.open(encoding="utf-8", newline="") as handle:
                mode_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["mode"] for row in mode_rows], ["taxa", "group", "node"])
            by_mode = {row["mode"]: row for row in mode_rows}
            self.assertEqual(by_mode["taxa"]["access_executed"], "YES")
            self.assertEqual(by_mode["taxa"]["collapse_executed"], "NO")
            self.assertEqual(by_mode["taxa"]["terminal_stage_used_for_final"], "access")
            self.assertEqual(by_mode["group"]["access_executed"], "NO")
            self.assertEqual(by_mode["group"]["collapse_executed"], "YES")
            self.assertEqual(by_mode["group"]["terminal_stage_used_for_final"], "collapse")
            self.assertEqual(by_mode["node"]["terminal_stage_used_for_final"], "generation")
            with final_merge_report.open(encoding="utf-8", newline="") as handle:
                final_merge_row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                final_merge_row["unique_sequences_in_merged_inputs"],
                final_merge_row["records_in_final_probe_set_fasta"],
            )
            self.assertEqual(final_merge_row["final_dedup_executed"], "YES")

            taxonomy_report = report_dir / "input_taxonomy_summary.tsv"
            self.assertEqual(
                input_taxonomy_statistics.main(
                    [
                        "--output-root",
                        str(output_root),
                        "--report",
                        str(taxonomy_report),
                    ]
                ),
                0,
            )
            with taxonomy_report.open(encoding="utf-8", newline="") as handle:
                taxonomy_row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(taxonomy_row["raw_input_genomes"], "13")
            self.assertEqual(taxonomy_row["excluded_by_length_qc"], "1")
            self.assertEqual(taxonomy_row["retained_genomes"], "12")
            self.assertEqual(taxonomy_row["raw_species"], "6")
            self.assertEqual(taxonomy_row["raw_genera"], "5")
            self.assertEqual(taxonomy_row["raw_families"], "2")
            self.assertEqual(taxonomy_row["retained_species"], "6")
            self.assertEqual(taxonomy_row["retained_genera"], "5")
            self.assertEqual(taxonomy_row["retained_families"], "2")

            access_report_dir = report_dir / "taxa_access_filters"
            self.assertEqual(
                access_reporting.main(
                    [
                        "--output-root",
                        str(output_root),
                        "--mode",
                        "taxa",
                        "--report-dir",
                        str(access_report_dir),
                    ]
                ),
                0,
            )
            with (access_report_dir / "taxa_access_filter_summary.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                access_summary = {
                    row["metric"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(
                access_summary["recomputed_assessed"],
                access_summary["actual_assessed"],
            )


if __name__ == "__main__":
    unittest.main()
