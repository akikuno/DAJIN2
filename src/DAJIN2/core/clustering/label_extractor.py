from __future__ import annotations

import uuid
from itertools import groupby
from pathlib import Path

from DAJIN2.core.clustering.clustering import return_labels
from DAJIN2.core.clustering.label_updator import relabel_with_consective_order
from DAJIN2.core.clustering.score_handler import annotate_score, make_score
from DAJIN2.core.clustering.strand_bias_handler import is_strand_bias
from DAJIN2.utils import io


def extract_labels(classif_sample, TEMPDIR, SAMPLE_NAME, CONTROL_NAME) -> list[dict[str]]:
    labels_result = []
    max_label = 1

    strand_bias = is_strand_bias(Path(TEMPDIR, CONTROL_NAME, "midsv", "control", f"{CONTROL_NAME}.jsonl"))
    classif_sample.sort(key=lambda x: x["ALLELE"])
    min_cluster_size = max(5, int(len(classif_sample) * 0.5 // 100))  # 0.5% of the samples

    for allele, group in groupby(classif_sample, key=lambda x: x["ALLELE"]):
        # Cache data to temporary files

        # For Insertion/Inversion allele
        path_control = Path(TEMPDIR, CONTROL_NAME, "midsv", allele, f"{SAMPLE_NAME}.jsonl")
        if not path_control.exists():
            path_control = Path(TEMPDIR, CONTROL_NAME, "midsv", allele, f"{CONTROL_NAME}.jsonl")

        unique_id = str(uuid.uuid4())
        path_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"tmp_{allele}_{unique_id}.jsonl")
        io.write_jsonl(data=group, file_path=path_sample)

        # Load mutation_loci and knockin_loci
        path_mutation_loci = Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", allele, "mutation_loci.pickle")
        path_knockin_loci = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", allele, "knockin.pickle")

        mutation_loci: list[set[str]] = io.load_pickle(path_mutation_loci)
        knockin_loci: set[int] = io.load_pickle(path_knockin_loci) if path_knockin_loci.exists() else set()

        # Skip clustering when the number of reads is too small or there is no mutation
        read_numbers = io.count_newlines(path_sample)
        is_no_mutation = all(m == set() for m in mutation_loci)
        if read_numbers < max(5, min_cluster_size) or is_no_mutation:
            labels_result.extend([max_label] * read_numbers)
            max_label += 1
            continue

        # Calculate scores to temporary files
        mutation_score: list[dict[str, float]] = make_score(path_sample, path_control, mutation_loci, knockin_loci)

        scores_sample = annotate_score(path_sample, mutation_score, mutation_loci)
        scores_control = annotate_score(path_control, mutation_score, mutation_loci, is_control=True)

        path_score_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"tmp_{allele}_score_sample_{unique_id}.jsonl")
        path_score_control = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"tmp_{allele}_score_control_{unique_id}.jsonl")
        io.write_jsonl(data=scores_sample, file_path=path_score_sample)
        io.write_jsonl(data=scores_control, file_path=path_score_control)

        # Extract labels
        labels = return_labels(path_score_sample, path_score_control, path_sample, strand_bias)
        labels_result += relabel_with_consective_order(labels, start=max_label)

        # Remove temporary files
        path_sample.unlink()
        path_score_sample.unlink()
        path_score_control.unlink()

        max_label = max(labels_result) + 1

    return labels_result
