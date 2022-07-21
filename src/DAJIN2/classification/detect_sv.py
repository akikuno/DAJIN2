def is_sv(CSSPLIT: str, kmer_size=100, window_size=50, threshold=50) -> bool:
    """Identify Structural variants (SVs).
    Structural variants are large genomic alterations, where large is typically (and somewhat arbitrarily) defined as encompassing at least 50 bp.

    Args:
        CSSPLIT (str): CSSPLIT from midsv conversion
        kmer_size (int, optional): Nucleotide sequence of a certain length. Defaults to 100.
        window_size (int, optional): Size of sliding window. Defaults to 50.
        threshold (int, optional): Threshold for the nucreotide number to be considered as structural variants. Defaults to 50.

    Raises:
        ValueError: kmer_size must be larger than window_size

    Returns:
        bool: Return True when a read is SV
    """
    if kmer_size < window_size:
        raise ValueError("kmer_size must be larger than window_size")
    test_list = CSSPLIT.split(",")
    is_SV = False
    for start in range(0, len(test_list) - window_size, window_size):
        end = start + kmer_size
        test_sliced = test_list[start:end]
        count_mutation = []
        for mutation in ["*", "-", "+", "N"]:
            count_mutation.append(sum(n.count(mutation) for n in test_sliced))
        count_mutation.append(sum(n.islower() for n in test_sliced))
        if any(c >= threshold for c in count_mutation):
            is_SV = True
        if is_SV:
            break
    return is_SV

