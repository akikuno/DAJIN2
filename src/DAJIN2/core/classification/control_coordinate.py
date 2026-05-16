from __future__ import annotations

import gzip
from collections import Counter, defaultdict
from pathlib import Path

from rapidfuzz.distance import Levenshtein

from DAJIN2.core.preprocess.alignment import mapping
from DAJIN2.core.preprocess.alignment.midsv_caller import (
    append_best_preset,
    convert_consecutive_indels,
    convert_flag_to_strand,
    filter_samples_by_n_proportion,
    replace_internal_n_to_d,
    transform_to_midsv_format,
)
from DAJIN2.core.preprocess.error_correction.sequence_error_fastq import (
    extract_qname_from_fastq_identifier,
    format_fastq_record,
)
from DAJIN2.utils import fileio, midsv_handler
from DAJIN2.utils.allele_handler import to_allele_key
from DAJIN2.utils.fastx_handler import read_fastq

PRESETS = ["map-ont", "map-hifi", "sr"]


def reference_signature(sequence: str) -> list[str]:
    return [f"={base}" for base in sequence]


def _insertion_payload(tag: str) -> str:
    return "".join(part[1:] for part in tag.split("|") if part.startswith("+"))


def _token_distance(read_token: str, allele_token: str) -> int:
    if read_token == allele_token:
        return 0
    if midsv_handler.is_n_tag(read_token) or midsv_handler.is_n_tag(allele_token):
        return 0
    if read_token.startswith("+") and allele_token.startswith("+"):
        return 1 + Levenshtein.distance(_insertion_payload(read_token), _insertion_payload(allele_token))
    if read_token[0] == allele_token[0]:
        return 1
    return 2


def score_midsv_against_signature(read_midsv: str | list[str], allele_signature: list[str]) -> int:
    read_tokens = read_midsv.split(",") if isinstance(read_midsv, str) else read_midsv
    shared_length = min(len(read_tokens), len(allele_signature))
    penalty = abs(len(read_tokens) - len(allele_signature)) * 2
    penalty += sum(_token_distance(read_tokens[i], allele_signature[i]) for i in range(shared_length))
    return -penalty


def signature_to_sequence(signature: list[str]) -> str:
    return midsv_handler.convert_midsvs_to_sequence(signature)


def _rank_alleles_by_sequence(
    read_midsv: str, allele_sequences: dict[str, str], allele_signatures: dict[str, list[str]]
) -> list[tuple[str, int]]:
    read_sequence = midsv_handler.convert_midsvs_to_sequence(read_midsv.split(","))
    scores = []
    for allele, allele_sequence in allele_sequences.items():
        sequence_score = -Levenshtein.distance(read_sequence, allele_sequence)
        token_score = score_midsv_against_signature(read_midsv, allele_signatures[allele])
        scores.append((allele, sequence_score, token_score))

    scores.sort(key=lambda item: (item[1], item[2], item[0]), reverse=True)
    return [(allele, sequence_score) for allele, sequence_score, _token_score in scores]


def rank_alleles(read_midsv: str, allele_signatures: dict[str, list[str]]) -> list[tuple[str, int]]:
    allele_sequences = {allele: signature_to_sequence(signature) for allele, signature in allele_signatures.items()}
    return _rank_alleles_by_sequence(read_midsv, allele_sequences, allele_signatures)


def classify_records_by_control_coordinate(
    midsv_records: list[dict], allele_signatures: dict[str, list[str]], no_filter: bool = False
) -> tuple[list[dict], dict[str, str]]:
    ranked_by_qname: dict[str, list[tuple[str, int]]] = {}
    record_by_qname: dict[str, dict] = {}
    allele_sequences = {allele: signature_to_sequence(signature) for allele, signature in allele_signatures.items()}

    for record in midsv_records:
        qname = record["QNAME"]
        ranked_by_qname[qname] = _rank_alleles_by_sequence(record["MIDSV"], allele_sequences, allele_signatures)
        record_by_qname[qname] = record

    assignments = _assign_alleles(ranked_by_qname, no_filter=no_filter)
    classified = []
    for qname, allele in assignments.items():
        record = record_by_qname[qname].copy()
        record["ALLELE"] = allele
        record["SCORE"] = ranked_by_qname[qname][0][1]
        classified.append(record)

    return classified, assignments


def _assign_alleles(ranked_by_qname: dict[str, list[tuple[str, int]]], no_filter: bool = False) -> dict[str, str]:
    primary_assignments = {qname: ranked[0][0] for qname, ranked in ranked_by_qname.items() if ranked}
    if no_filter:
        return primary_assignments

    allele_counts = Counter(primary_assignments.values())
    minor_alleles = {allele for allele, count in allele_counts.items() if count < 5}
    if not minor_alleles:
        return primary_assignments

    major_alleles = set(allele_counts) - minor_alleles
    if not major_alleles:
        return primary_assignments

    assignments = primary_assignments.copy()
    most_common_major = max(major_alleles, key=lambda allele: allele_counts[allele])
    for qname, allele in primary_assignments.items():
        if allele not in minor_alleles:
            continue
        assignments[qname] = next(
            (candidate for candidate, _score in ranked_by_qname[qname] if candidate in major_alleles),
            most_common_major,
        )
    return assignments


def load_control_midsv_records(tempdir: Path, sample_name: str) -> list[dict]:
    control_key = to_allele_key("control")
    path_midsv = Path(tempdir, sample_name, "midsv", control_key, f"{sample_name}_midsv.jsonl")
    return list(fileio.read_jsonl(path_midsv))


def build_allele_signatures(ARGS) -> dict[str, list[str]]:
    path_cache = Path(ARGS.tempdir, "cache", "control_coordinate", "allele_signatures.pickle")
    path_cache.parent.mkdir(parents=True, exist_ok=True)
    if path_cache.exists():
        return fileio.load_pickle(path_cache)

    signatures = {}
    control_sequence = ARGS.fasta_alleles["control"]
    for allele, sequence in ARGS.fasta_alleles.items():
        path_consensus_midsv = Path(
            ARGS.tempdir, ARGS.sample_name, "midsv", f"consensus_{to_allele_key(allele)}_midsv.jsonl"
        )
        if allele == "control":
            signatures[allele] = reference_signature(control_sequence)
        elif path_consensus_midsv.exists():
            signatures[allele] = list(fileio.read_jsonl(path_consensus_midsv))
        else:
            signatures[allele] = align_allele_to_control_signature(ARGS, allele, sequence)

    fileio.save_pickle(signatures, path_cache)
    return signatures


def align_allele_to_control_signature(ARGS, allele: str, sequence: str) -> list[str]:
    control_key = to_allele_key("control")
    allele_key = to_allele_key(allele)
    path_reference = Path(ARGS.tempdir, ARGS.control_name, "fasta", f"{control_key}.fasta")
    path_query = Path(ARGS.tempdir, "cache", "control_coordinate", "allele_queries", f"{allele_key}.fasta")
    path_sam = Path(ARGS.tempdir, "cache", "control_coordinate", "allele_sam", f"{allele_key}.sam")
    path_query.parent.mkdir(parents=True, exist_ok=True)
    path_sam.parent.mkdir(parents=True, exist_ok=True)
    path_query.write_text(f">{allele}\n{sequence}\n")

    sam_lines = list(mapping.to_sam(path_reference, path_query, preset="map-ont", threads=ARGS.threads))
    path_sam.write_text("\n".join(sam_lines))
    if sum(1 for line in sam_lines if not line.startswith("@")) == 0:
        return reference_signature(ARGS.fasta_alleles["control"])

    midsv_records = transform_to_midsv_format(path_sam)
    midsv_records = replace_internal_n_to_d(midsv_records, ARGS.fasta_alleles["control"])
    midsv_records = convert_flag_to_strand(midsv_records)
    midsv_records = filter_samples_by_n_proportion(midsv_records)
    midsv_records = convert_consecutive_indels(midsv_records)
    midsv_records = list(append_best_preset(midsv_records, {allele: "map-ont"}))
    if not midsv_records:
        return reference_signature(ARGS.fasta_alleles["control"])
    return midsv_records[0]["MIDSV"].split(",")


def split_fastq_by_assignment(path_fastq: Path, assignments: dict[str, str], path_output_dir: Path) -> dict[str, Path]:
    path_output_dir.mkdir(parents=True, exist_ok=True)
    reads_by_allele = defaultdict(list)
    for read in read_fastq(path_fastq):
        qname = extract_qname_from_fastq_identifier(read["identifier"])
        allele = assignments.get(qname)
        if allele is None:
            continue
        reads_by_allele[allele].append(read)

    paths_by_allele = {}
    for allele, reads in reads_by_allele.items():
        allele_key = to_allele_key(allele)
        path_output = Path(path_output_dir, f"{allele_key}.fastq.gz")
        with gzip.open(path_output, "wt") as f:
            for read in reads:
                f.write(format_fastq_record(read))
        paths_by_allele[allele] = path_output
    return paths_by_allele


def write_assigned_sam_files(
    ARGS,
    paths_fastq_by_allele: dict[str, Path],
    paths_fasta_by_allele: dict[str, Path],
    name: str,
    is_control_sv: bool = False,
) -> None:
    for allele, path_fastq in paths_fastq_by_allele.items():
        path_fasta = paths_fasta_by_allele[allele]
        allele_key = to_allele_key(allele)
        path_sam_directory = Path(ARGS.tempdir, name, "sam", allele_key)
        path_sam_directory.mkdir(parents=True, exist_ok=True)
        for preset in PRESETS:
            sam = mapping.to_sam(path_fasta, path_fastq, preset=preset, threads=ARGS.threads)
            if is_control_sv:
                path_sam_file = Path(path_sam_directory, f"{ARGS.sample_name}_{preset}.sam")
            else:
                path_sam_file = Path(path_sam_directory, f"{preset}.sam")
            path_sam_file.write_text("\n".join(sam))


def load_classified_midsv_for_assignments(ARGS, assignments: dict[str, str]) -> list[dict]:
    qnames_by_allele = defaultdict(set)
    for qname, allele in assignments.items():
        qnames_by_allele[allele].add(qname)

    classified = []
    for allele, qnames in qnames_by_allele.items():
        allele_key = to_allele_key(allele)
        path_midsv = Path(ARGS.tempdir, ARGS.sample_name, "midsv", allele_key, f"{ARGS.sample_name}_midsv.jsonl")
        if not path_midsv.exists():
            continue
        for record in fileio.read_jsonl(path_midsv):
            if record["QNAME"] not in qnames:
                continue
            record["ALLELE"] = allele
            classified.append(record)
    return classified
