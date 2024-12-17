from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path

import cstag
import numpy as np
from Bio.Align import PairwiseAligner

from DAJIN2.core.clustering.clustering import optimize_labels
from DAJIN2.core.preprocess.sv_handler import add_unique_allele_keys, extract_unique_sv, save_cstag, save_fasta
from DAJIN2.utils import io
from DAJIN2.utils.cssplits_handler import convert_cssplits_to_cstag

###############################################################################
# Detect sequence errors
###############################################################################


def convert_sv_tag_insertion(cssplits: list[str]) -> str:
    """Create a tag with I for bases that reads with insertion and M for bases that are mapped"""
    im_tags = []
    for tag in cssplits.split(","):
        if tag.startswith("+") and tag.count("+") > 10:
            ins = "I" * tag.count("|")
            im_tags.append(ins + "M")
        else:
            im_tags.append("M")

    return "".join(im_tags)


# TODO: def convert_sv_tag_deletion(cssplits: list[str]) -> str:
# TODO: def convert_sv_tag_inversion(cssplits: list[str]) -> str:


def convert_sv_tag(cssplits: list[str], sv_type: str) -> str:
    if sv_type == "insertion":
        return convert_sv_tag_insertion(cssplits)
    # TODO


###############################################################################
# extract_sv_features
###############################################################################


def _find_sv_start_index(sv_tag: str, sv_type: str) -> list[int]:
    if sv_type == "insertion":
        key = "I"
    elif sv_type == "deletion":
        key = "D"
    elif sv_type == "inversion":
        key = "V"

    pattern = re.compile(rf"{key}+")
    return [match.start() for match in pattern.finditer(sv_tag)]


def _extract_sv_start_index(sv_tags: list[str], sv_type: str) -> list[int]:
    sv_start_index = []
    for tag in sv_tags:
        for i in _find_sv_start_index(tag, sv_type):
            sv_start_index.append(i)
    return sorted(sv_start_index)


def _group_sv_start_index(sv_tags: list[list[str]], sv_type: str) -> list[list[int]]:
    sv_start_index = _extract_sv_start_index(sv_tags, sv_type)

    groups = []
    current_group = []
    for key in sv_start_index:
        if not current_group or key - current_group[-1] <= 5:
            current_group.append(key)
        else:
            groups.append(current_group)
            current_group = [key]

    if current_group:
        groups.append(current_group)

    return [group for group in groups if len(group) > 10]


def define_index_converter(sv_tags: list[list[str]], sv_type: str) -> dict[int, int]:
    index_group = _group_sv_start_index(sv_tags, sv_type)
    index_converter = {}
    for group in index_group:
        value = group[0]
        index_converter |= {g: value for g in group}
    return index_converter


def annotate_merged_sv_start_index(
    sv_tags: list[str], sv_type: str, index_converter: dict[int, int]
) -> list[list[int]]:
    sv_start_index = []
    for tag in sv_tags:
        start_index = []
        for i in _find_sv_start_index(tag, sv_type):
            start_index.append(index_converter.get(i, None))
        sv_start_index.append(start_index)
    return sv_start_index


def extract_sv_start_index(sv_tags_index: list[list[int]], sv_index: list[int]) -> list[list[int]]:
    X_index = []
    for index in sv_tags_index:
        x_index = [0] * len(sv_index)
        for i in index:
            if i:
                x_index[sv_index.index(i)] = 1
        X_index.append(x_index)
    return X_index


def _find_sv_start_index_and_size(sv_tag: str, sv_type: str) -> list[int]:
    if sv_type == "insertion":
        key = "I"
    elif sv_type == "deletion":
        key = "D"
    elif sv_type == "inversion":
        key = "V"

    pattern = re.compile(rf"{key}+")
    return [[match.start(), match.end() - match.start()] for match in pattern.finditer(sv_tag)]


def annotate_size_by_merged_sv_start_index(
    sv_tags: list[str], sv_type: str, index_converter: dict[int, int]
) -> list[list[int]]:
    sv_start_index_size = []
    for tag in sv_tags:
        start_index_size = []
        for i, size in _find_sv_start_index_and_size(tag, sv_type):
            start_index_size.append([index_converter.get(i, None), size])
        sv_start_index_size.append(start_index_size)
    return sv_start_index_size


def extract_sv_size(sv_tags_size: list[list[list[int]]], sv_index: list[int]) -> list[list[int]]:
    X_size = []
    for index_size in sv_tags_size:
        x_size = [0] * len(sv_index)
        for i, size in index_size:
            if i:
                x_size[sv_index.index(i)] = size
        X_size.append(x_size)
    return X_size


def extract_sv_features(sv_tags: list[str], sv_type: str, sv_index: list[int]) -> np.ndarray[np.float64]:
    index_converter = define_index_converter(sv_tags, sv_type)
    sv_tags_index = annotate_merged_sv_start_index(sv_tags, sv_type, index_converter)
    sv_tags_size = annotate_size_by_merged_sv_start_index(sv_tags, sv_type, index_converter)

    X_index = extract_sv_start_index(sv_tags_index, sv_index)
    X_size = extract_sv_size(sv_tags_size, sv_index)

    return np.concatenate([X_index, X_size], axis=1)


###############################################################################
# ToDO 変異部についてコンセンサス配列を入手する
###############################################################################


def _group_cssplits(cssplits: list[str]) -> Iterator[str]:
    current_group = [cssplits[0]]
    current_start_index = 0

    for i in range(1, len(cssplits)):
        # 現在の要素のprefixと前の要素のprefixが同じ場合
        if cssplits[i][0] == cssplits[i - 1][0]:
            current_group.append(cssplits[i])
        else:
            yield (current_start_index, ",".join(current_group))  # グループを結果に追加
            current_group = [cssplits[i]]  # 新しいグループを開始
            current_start_index = i

    yield (i, ",".join(current_group))  # 最後のグループを追加


def _filter_valid_cssplits(cssplits, index_converter):
    """有効な cssplits をフィルタリングする"""
    # TODO: Deletion, Insertion, Inversionで適切な処理を変更する
    return [
        (i, cssplit)
        for i, cssplit in cssplits
        if cssplit.startswith("+") and cssplit.count("+") > 10 and i in index_converter
    ]


def _convert_insertion_cstag_to_cssplit(cs_tag: str) -> str:
    # TODO: Deletion, Insertion, Inversionで適切な処理を変更する
    pattern = re.compile(r"([\+\-=])([a-zA-Z]+)(.*)")
    match = pattern.match(cs_tag)
    if not match:
        return cs_tag  # マッチしない場合はそのまま返す
    prefix, letters, suffix = match.groups()  # グループを取得
    # 文字ごとにprefixを付けて'|'で結合
    transformed = "|".join(f"{prefix}{char}" for char in letters)

    # 最後の部分を追加
    return f"{transformed}{suffix}".upper()


# すべての配列をシードとアラインメントし、その結果を保存
def _align_sequences(seed: str, sequences: str):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligned_sequences = []
    for seq in sequences:
        alignment = aligner.align(seed, seq)[0]  # 最初のアラインメント結果を取得
        aligned_seq = []
        for a, b in zip(alignment.target, alignment.query):
            if a == "-":
                aligned_seq.append("-")
            else:
                aligned_seq.append(b)
        aligned_sequences.append("".join(aligned_seq))
    return aligned_sequences


# アラインメント結果からコンセンサス配列を生成
def _build_consensus(aligned_seqs):
    max_len = max(len(seq) for seq in aligned_seqs)
    padded_seqs = [seq.ljust(max_len, "-") for seq in aligned_seqs]  # ギャップで揃える
    aligned_matrix = np.array([list(seq) for seq in padded_seqs])
    consensus = []
    for col in aligned_matrix.T:  # 列ごとに塩基をチェック
        unique, counts = np.unique(col, return_counts=True)
        most_common = unique[np.argmax(counts)]  # 最も頻度の高い塩基を選定
        consensus.append(most_common if most_common != "-" else "-")
    return "".join(consensus).replace("-", "")  # 最終的にギャップを取り除く


def _call_consensus(sequences: list[str]) -> str:
    # 最初に最も長いシード配列を選定
    seed = max(sequences, key=len)
    aligned_sequences = _align_sequences(seed, sequences)
    return _build_consensus(aligned_sequences)


def _generate_consensus_for_label(cssplits_subset, index_converter, sv_type: str) -> dict[int, str]:
    """cssplits_subset からインデックスごとのコンセンサスを生成する"""
    csv_group = defaultdict(list)

    for cssplits in cssplits_subset:
        for i, cssplit in _group_cssplits(cssplits):
            if_valid = _filter_valid_cssplits([(i, cssplit)], index_converter)  # TODO
            for i, valid_cssplit in if_valid:
                cs_tags = convert_cssplits_to_cstag([valid_cssplit])
                csv_group[index_converter[i]].append(cs_tags)

    index_consensus = {}
    for key, value in csv_group.items():
        if sv_type == "insertion":
            ins_seqs, end_seqs = zip(*[cstag.split(tag) for tag in value])
            consensus_ins = _call_consensus(ins_seqs) + "|" + _call_consensus(end_seqs)
            index_consensus[key] = _convert_insertion_cstag_to_cssplit(consensus_ins)  # TODO
        else:
            index_consensus[key] = _call_consensus(value)
    return index_consensus


def _replace_cssplits_with_consensus(cssplits_subset, index_consensus, index_converter):
    """cssplits_subset の中でコンセンサスに置き換える"""
    for cssplits in cssplits_subset:
        for i, cssplit in _group_cssplits(cssplits):
            if_valid = _filter_valid_cssplits([(i, cssplit)], index_converter)
            for i, _ in if_valid:
                cssplits[i] = index_consensus[index_converter[i]]


###############################################################################
# cstagとfastaをラベルごとに生成する
###############################################################################


def generate_cstag_and_fasta(
    path_midsv: Path, labels: list[int], index_converter: dict[int, int], sv_type: str
) -> tuple[dict[int, str], dict[int, str]]:
    """cstag と fasta をラベルごとに生成する"""

    midsv_sample = io.read_jsonl(path_midsv)
    cssplits_iter: Iterator[list[list[str]]] = (m["CSSPLIT"].split(",") for m in midsv_sample)

    cssplits_by_label = defaultdict(list)
    for label, cssplits in zip(labels, cssplits_iter):
        cssplits_by_label[label].append(cssplits)

    cstag_by_label = {}
    fasta_by_label = {}

    for label, cssplits in cssplits_by_label.items():
        cssplits_subset = cssplits[:100]

        # コンセンサスの生成
        index_consensus = _generate_consensus_for_label(cssplits_subset, index_converter, sv_type)
        _replace_cssplits_with_consensus(cssplits_subset, index_consensus, index_converter)

        # コンセンサス cstag と fasta の生成
        cstag_subset = [convert_cssplits_to_cstag(tag) for tag in cssplits_subset]
        positions = [1] * len(cstag_subset)
        cstag_consensus = cstag.consensus(cstag_subset, positions)

        cstag_by_label[label] = cstag_consensus
        fasta_by_label[label] = cstag.to_sequence(cstag_consensus)

    return cstag_by_label, fasta_by_label


###############################################################################
# main
###############################################################################


def detect_sv_alleles(TEMPDIR: Path, SAMPLE_NAME: str, CONTROL_NAME: str, FASTA_ALLELES: dict, sv_type: str) -> None:
    path_midsv_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control", f"{SAMPLE_NAME}.jsonl")
    midsv_sample = io.read_jsonl(path_midsv_sample)
    sv_tags_sample = [convert_sv_tag(m["CSSPLIT"], sv_type) for m in midsv_sample]
    MUTATION_LOCI = io.load_pickle(
        Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", "control", "mutation_loci.pickle")
    )  # TODO: MUTATION_LOCIをコンセンサスに反映させる

    path_midsv_control = Path(TEMPDIR, CONTROL_NAME, "midsv", "control", f"{CONTROL_NAME}.jsonl")
    midsv_control = io.read_jsonl(path_midsv_control)
    sv_tags_control = [convert_sv_tag(m["CSSPLIT"], sv_type) for m in midsv_control][:1000]

    index_converter = define_index_converter(sv_tags_sample, sv_type)
    sv_index = sorted(set(index_converter.values()))

    sv_features_sample = extract_sv_features(sv_tags_sample, sv_type, sv_index)
    sv_features_control = extract_sv_features(sv_tags_control, sv_type, sv_index)

    #######################################################
    # Clustering
    #######################################################
    X = np.concatenate([sv_features_sample, sv_features_control])
    coverage_control = len(sv_features_control)
    coverage_sample = len(sv_features_sample)

    labels = optimize_labels(X, coverage_sample, coverage_control)

    #######################################################
    # Generate cstag consensus
    #######################################################

    cstag_by_label, fasta_by_label = generate_cstag_and_fasta(path_midsv_sample, labels, index_converter, sv_type)

    #######################################################
    # Remove similar alleles to user's alleles, or clustered alleles
    #######################################################
    fasta_by_label = extract_unique_sv(fasta_by_label, FASTA_ALLELES)
    cstag_by_label = {label: cs_tag for label, cs_tag in cstag_by_label.items() if label in fasta_by_label}

    #######################################################
    # Output cstag and fasta
    #######################################################
    cstag_by_label = add_unique_allele_keys(cstag_by_label, FASTA_ALLELES, key=sv_type)
    fasta_by_label = add_unique_allele_keys(fasta_by_label, FASTA_ALLELES, key=sv_type)

    save_cstag(TEMPDIR, SAMPLE_NAME, cstag_by_label)
    save_fasta(TEMPDIR, SAMPLE_NAME, fasta_by_label)
