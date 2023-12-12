from __future__ import annotations

from DAJIN2.utils import io, config

config.set_warnings()

from pathlib import Path
import random
from itertools import groupby

from DAJIN2.core.clustering.score_handler import make_score, annotate_score
from DAJIN2.core.clustering.label_updator import relabel_with_consective_order
from DAJIN2.core.clustering.strand_bias_handler import is_strand_bias
from DAJIN2.core.clustering.clustering import return_labels


def extract_labels(classif_sample, TEMPDIR, SAMPLE_NAME, CONTROL_NAME) -> list[dict[str]]:
    labels_all = []
    max_label = 0
    strand_bias = is_strand_bias(Path(TEMPDIR, CONTROL_NAME, "midsv", "control.json"))
    classif_sample.sort(key=lambda x: x["ALLELE"])
    for allele, group in groupby(classif_sample, key=lambda x: x["ALLELE"]):
        path_mutation_loci = Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", f"{allele}.pickle")
        mutation_loci: list[set[str]] = io.load_pickle(path_mutation_loci)
        if all(m == set() for m in mutation_loci):
            max_label += 1
            labels_all.extend([max_label] * len(classif_sample))
            continue

        path_knockin_loci = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", f"{allele}.pickle")
        knockin_loci: set[int] = io.load_pickle(path_knockin_loci) if path_knockin_loci.exists() else set()

        RANDOM_INT = random.randint(0, 10**10)
        path_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"{allele}_{RANDOM_INT}.json")
        if Path(TEMPDIR, CONTROL_NAME, "midsv", f"{allele}.json").exists():
            path_control = Path(TEMPDIR, CONTROL_NAME, "midsv", f"{allele}.json")
        else:
            path_control = Path(TEMPDIR, CONTROL_NAME, "midsv", f"{allele}_{SAMPLE_NAME}.json")
        io.write_jsonl(data=group, file_path=path_sample)

        """Prepare and write clustering data to temporary files."""
        mutation_score: list[dict[str, float]] = make_score(path_sample, path_control, mutation_loci, knockin_loci)

        scores_sample = annotate_score(path_sample, mutation_score, mutation_loci)
        scores_control = annotate_score(path_control, mutation_score, mutation_loci, is_control=True)

        path_score_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"{allele}_score_{RANDOM_INT}.json")
        path_score_control = Path(TEMPDIR, CONTROL_NAME, "clustering", f"{allele}_score_{RANDOM_INT}.json")
        io.write_jsonl(data=scores_sample, file_path=path_score_sample)
        io.write_jsonl(data=scores_control, file_path=path_score_control)

        """Extract labels."""
        labels = return_labels(path_score_sample, path_score_control, path_sample, strand_bias)
        labels_reordered = relabel_with_consective_order(labels, start=max_label)

        max_label = max(labels_reordered)
        labels_all.extend(labels_reordered)

        """Remove temporary files."""
        path_sample.unlink()
        path_score_sample.unlink()
        path_score_control.unlink()

    return labels_all
