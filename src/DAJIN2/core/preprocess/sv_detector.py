from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Iterator
from itertools import chain
from pathlib import Path

import cstag
import numpy as np

from DAJIN2.core.clustering.clustering import optimize_labels
from DAJIN2.core.preprocess.sv_handler import add_unique_allele_keys, extract_unique_sv, save_cstag, save_fasta
from DAJIN2.utils import io
from DAJIN2.utils.cssplits_handler import convert_cssplits_to_cstag

###############################################################################
# Clusteringにつかう特徴量を抽出する
# 特徴量：SVの開始のIndexとSVの大きさ
###############################################################################


def get_index_of_sv(misdv_string: str, sv_type: str, mutation_loci: list[set[str]]) -> list[int]:
    """Return a list with the start index of SV"""

    items = misdv_string.split(",")
    index_of_sv = []
    previous_midsv = items[0]

    for i, item in enumerate(items):
        if sv_type in ("deletion", "inversion"):
            if sv_type == "deletion":
                is_sv_item = item.startswith("-") and "-" in mutation_loci[i]
                is_prev_sv_item = not previous_midsv.startswith("-")
            else:  # inversion
                is_sv_item = item.islower()
                is_prev_sv_item = not previous_midsv.islower()

            if is_sv_item and is_prev_sv_item:
                index_of_sv.append(i)

        elif sv_type == "insertion":
            if item.startswith("+") and "+" in mutation_loci[i]:
                index_of_sv.append(i)

        previous_midsv = item
    return index_of_sv


def group_similar_indices(index_of_sv: list[int], distance: int = 5, min_coverage: float = 5.0) -> list[list[int]]:
    """Return a list of groups of similar indices"""
    index_of_sv.sort()
    groups = []
    current_group = []
    for key in index_of_sv:
        if not current_group or key - current_group[-1] <= distance:
            current_group.append(key)
        else:
            groups.append(current_group)
            current_group = [key]

    if current_group:
        groups.append(current_group)

    return [group for group in groups if len(group) > min_coverage]


def define_index_converter(index_grouped: list[list[int]]) -> dict[int, int]:
    index_converter = {}
    for group in index_grouped:
        value = Counter(group).most_common()[0][0]
        index_converter |= {g: value for g in group}
    return index_converter


def get_index_and_sv_size(
    misdv_string: str, sv_type: str, mutation_loci: list[set[str]], index_converter: dict[int, int]
) -> dict[int, int]:
    """Return a dictionary with the start index and SV size"""

    items = misdv_string.split(",")
    result = {}
    sv_size = 0
    start_index = None
    previous_midsv = items[0]

    for i, item in enumerate(items):
        if sv_type in ("deletion", "inversion"):
            if sv_type == "deletion":
                is_sv_item = item.startswith("-") and "-" in mutation_loci[i]
                is_prev_sv_item = previous_midsv.startswith("-")
            else:  # inversion
                is_sv_item = item.islower()
                is_prev_sv_item = previous_midsv.islower()

            # Get the start index of SV
            if is_sv_item and not is_prev_sv_item and i in index_converter:
                start_index = index_converter[i]

            # Calculate the size of SV
            if is_sv_item:
                sv_size += 1
            elif start_index:
                result.update({start_index: sv_size})
                start_index = None
                sv_size = 0

        elif sv_type == "insertion":
            if item.startswith("+") and "+" in mutation_loci[i] and i in index_converter:
                result.update({index_converter[i]: item.count("|")})

        previous_midsv = item
    return result


def extract_features(index_and_sv_size: list[dict[int, int]], all_sv_index: set[str]) -> np.ndarray:
    """補正後のIndexにあるSVの大きさを特徴量として用いる"""
    sv_size_features = []

    for record in index_and_sv_size:
        sv_size = {i: 0 for i in all_sv_index}
        sv_size |= record
        sv_size_features.append(sv_size)

    index_of_sv = sorted({key for record in sv_size_features for key in record.keys()})
    return np.array([[record.get(i, 0) for i in index_of_sv] for record in sv_size_features])


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


def _filter_valid_cssplits(cssplits: list[tuple[int, str]], index_converter: dict[int, int]) -> list[tuple[int, str]]:
    """有効な cssplits をフィルタリングする"""
    # TODO: Deletion, Insertion, Inversionで適切な処理を変更する
    return [
        (i, cssplit)
        for i, cssplit in cssplits
        if cssplit.startswith("+") and cssplit.count("+") > 10 and i in index_converter
    ]


def _convert_insertion_cstag_to_cssplit(cs_tag: str) -> str:
    pattern = re.compile(r"([\+\-=])([a-zA-Z]+)(.*)")
    match = pattern.match(cs_tag)
    if not match:
        return cs_tag  # マッチしない場合はそのまま返す
    prefix, letters, suffix = match.groups()  # グループを取得
    # 文字ごとにprefixを付けて'|'で結合
    transformed = "|".join(f"{prefix}{char}" for char in letters)

    # 最後の部分を追加
    return f"{transformed}{suffix}".upper()


# ---------------------------------------------------------
# call_consensus
# ---------------------------------------------------------


# すべての配列をシードとアラインメントし、その結果を保存
def _align_sequences(seed: str, sequences: list[str]) -> list[str]:
    """
    Align all sequences to the seed sequence and return the aligned sequences.

    Args:
        seed (str): The seed sequence to align against.
        sequences (list[str]): The list of sequences to be aligned.

    Returns:
        list[str]: The list of aligned sequences.
    """
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
def _build_consensus(aligned_seqs: list[str]) -> str:
    """
    Generate a consensus sequence from the aligned sequences.

    Args:
        aligned_seqs (list[str]): The list of aligned sequences.

    Returns:
        str: The consensus sequence.
    """
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
            index_consensus[key] = _convert_insertion_cstag_to_cssplit(consensus_ins)
        else:
            index_consensus[key] = _call_consensus(value)

    return index_consensus


def _replace_cssplits_with_consensus(
    cssplits_subset: list[list[str]], index_consensus: dict[int, str], index_converter: dict[int, int]
):
    """cssplits_subset の中でコンセンサスに置き換える"""
    cssplits_replaced = cssplits_subset.copy()
    for cssplits in cssplits_replaced:
        for i, cssplit in _group_cssplits(cssplits):
            if_valid = _filter_valid_cssplits([(i, cssplit)], index_converter)
            for i, _ in if_valid:
                cssplits[i] = index_consensus[index_converter[i]]

    return cssplits_replaced


def _correct_sequence_error(
    cssplits_subset: list[list[str]], mutation_loci: list[set[str]], fasta_control: str
) -> list[list[str]]:
    cssplit_replaced = cssplits_subset.copy()

    for cssplits in cssplit_replaced:
        for i, cs in enumerate(cssplits):
            if cs[0] in mutation_loci[i]:
                continue
            cssplits[i] = "=" + fasta_control[i]

    return cssplit_replaced


###############################################################################
# cstagとfastaをラベルごとに生成する
###############################################################################


def generate_cstag_and_fasta(
    path_midsv: Path,
    mutation_loci: list[set[str]],
    fasta_control: str,
    labels: list[int],
    index_converter: dict[int, int],
    sv_type: str,
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
        cssplits_consensus = _replace_cssplits_with_consensus(cssplits_subset, index_consensus, index_converter)
        cssplits_consensus = _correct_sequence_error(cssplits_subset, mutation_loci, fasta_control)
        # コンセンサス cstag と fasta の生成
        cstag_subset = [convert_cssplits_to_cstag(tag) for tag in cssplits_consensus]
        positions = [1] * len(cstag_subset)
        cstag_consensus = cstag.consensus(cstag_subset, positions)

        cstag_by_label[label] = cstag_consensus
        fasta_by_label[label] = cstag.to_sequence(cstag_consensus)

    return cstag_by_label, fasta_by_label


###############################################################################
# main
###############################################################################


def detect_sv_alleles(TEMPDIR: Path, SAMPLE_NAME: str, CONTROL_NAME: str, FASTA_ALLELES: dict, sv_type: str) -> None:
    #######################################################
    # Load data
    #######################################################
    MUTATION_LOCI: list[set[str]] = io.load_pickle(
        Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", "control", "mutation_loci.pickle")
    )

    path_midsv_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control", f"{SAMPLE_NAME}.jsonl")
    path_midsv_control = Path(TEMPDIR, CONTROL_NAME, "midsv", "control", f"{CONTROL_NAME}.jsonl")

    coverage_sample = io.count_newlines(path_midsv_sample)
    coverage_control = io.count_newlines(path_midsv_control)

    #######################################################
    # Extract features: SV start index and SV size
    #######################################################
    index_of_sv_sample = [
        get_index_of_sv(m["CSSPLIT"], sv_type, MUTATION_LOCI) for m in io.read_jsonl(path_midsv_sample)
    ]
    index_of_sv_sample = sorted(chain.from_iterable(index_of_sv_sample))

    index_of_sv_control = [
        get_index_of_sv(m["CSSPLIT"], sv_type, MUTATION_LOCI) for m in io.read_jsonl(path_midsv_control)
    ][:1000]
    index_of_sv_control = sorted(chain.from_iterable(index_of_sv_control))

    index_grouped_sample = group_similar_indices(
        index_of_sv_sample, distance=5, min_coverage=max(5, coverage_sample * 0.05)
    )
    index_grouped_control = group_similar_indices(
        index_of_sv_control, distance=5, min_coverage=max(5, coverage_control * 0.05)
    )

    index_converter_sample = define_index_converter(index_grouped_sample)
    index_converter_control = define_index_converter(index_grouped_control)

    index_and_sv_size_sample = [
        get_index_and_sv_size(m["CSSPLIT"], sv_type, MUTATION_LOCI, index_converter_sample)
        for m in io.read_jsonl(path_midsv_sample)
    ]
    index_and_sv_size_control = [
        get_index_and_sv_size(m["CSSPLIT"], sv_type, MUTATION_LOCI, index_converter_control)
        for m in io.read_jsonl(path_midsv_control)
    ]

    all_sv_index = {
        key
        for record in index_and_sv_size_sample + index_and_sv_size_control
        if isinstance(record, dict)
        for key in record.keys()
    }

    X_sample = extract_features(index_and_sv_size_sample, all_sv_index)
    X_control = extract_features(index_and_sv_size_control, all_sv_index)

    #######################################################
    # Clustering
    #######################################################

    X = np.vstack([X_sample, X_control])
    labels = optimize_labels(X, coverage_sample, coverage_control)

    #######################################################
    # Call consensus
    #######################################################
    cstag_by_label = {}
    fasta_by_label = {}

    if sv_type == "deletion" or sv_type == "inversion":
        ###################################################
        # 各IndexにおけるSV sizeのコンセンサスを、コントロールのアレルに反映させる
        ###################################################

        sv_index_by_label = defaultdict(set)
        for label, record in zip(labels, index_and_sv_size_sample):
            sv_index_by_label[label] |= record.keys()

        index_and_sv_size_by_label = defaultdict(list)
        for label, record in zip(labels, index_and_sv_size_sample):
            sv_index = sv_index_by_label[label]
            sv_size = {i: 0 for i in sv_index}
            sv_size |= record
            index_and_sv_size_by_label[label].append(sv_size)

        ###################################################
        # SV sizeのコンセンサスを取得
        ###################################################
        consensus_index_and_sv_size_by_label = defaultdict(dict)
        for label, index_and_sv_size in index_and_sv_size_by_label.items():
            frequency_counts = defaultdict(lambda: defaultdict(int))

            # Count the frequencies
            for record in index_and_sv_size:
                for key, value in record.items():
                    frequency_counts[key][value] += 1

            # Extract the most frequent value for each key
            consensus_index_and_sv_size_by_label[label] = {
                key: max(values, key=values.get) for key, values in frequency_counts.items()
            }

        ###################################################
        # ControlのアレルにSVを反映させる
        ###################################################
        def modify_tags(midsv_tags: list[str], start: int, end: int, modify_func: callable) -> list[str]:
            return [tag if not (start <= i <= end) else modify_func(tag) for i, tag in enumerate(midsv_tags)]

        for label, consensus_index_and_sv_size in consensus_index_and_sv_size_by_label:
            midsv_tags_control = ["=" + s for s in list(FASTA_ALLELES["control"])]
            for start, sv_size in consensus_index_and_sv_size.items():
                if sv_size == 0:
                    continue
                end = start + sv_size
                if sv_type == "deletion":
                    midsv_tags_control = modify_tags(midsv_tags_control, start, end, lambda x: "-" + x[1:])
                elif sv_type == "inversion":
                    midsv_tags_control = modify_tags(midsv_tags_control, start, end, lambda x: x.lower())

            cstag_by_label[label] = convert_cssplits_to_cstag(midsv_tags_control)
            fasta_by_label[label] = cstag.to_sequence(cstag_by_label[label])

    else:  # sv_type == "insertion"
        cstag_by_label, fasta_by_label = generate_cstag_and_fasta(
            path_midsv_sample, MUTATION_LOCI, fasta_control, labels, index_converter_sample, sv_type
        )
    # TODO ---------------------

    for key, value in index_and_sv_size.items():
        if key not in consensus_index_and_sv_size_by_label[label]:
            consensus_index_and_sv_size_by_label[label][key] = value

    cssplits_iter: Iterator[list[list[str]]] = (m["CSSPLIT"].split(",") for m in io.read_jsonl(path_midsv_sample))

    cssplits_by_label = defaultdict(list)
    for label, cssplits in zip(labels, cssplits_iter):
        cssplits_by_label[label].append(cssplits)


def old_detect_sv_alleles(
    TEMPDIR: Path, SAMPLE_NAME: str, CONTROL_NAME: str, FASTA_ALLELES: dict, sv_type: str
) -> None:
    return
    path_midsv_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", "control", f"{SAMPLE_NAME}.jsonl")
    midsv_sample = io.read_jsonl(path_midsv_sample)
    sv_tags_sample = [convert_sv_tag(m["CSSPLIT"], sv_type) for m in midsv_sample]
    MUTATION_LOCI: list[set[str]] = io.load_pickle(
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
    cstag_by_label, fasta_by_label = generate_cstag_and_fasta(
        path_midsv_sample, MUTATION_LOCI, FASTA_ALLELES["control"], labels, index_converter, sv_type
    )

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
