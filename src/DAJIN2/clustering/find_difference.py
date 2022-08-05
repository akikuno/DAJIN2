from __future__ import annotations
import scipy.stats as st


def screen_different_loci(sample_cssplit: list[dict], control_cssplit: list[dict]) -> list[int]:
    """_summary_

    Args:
        path_sample_cssplit (Union[Path, str]): Path of sample's jsonl containing CSSPLIT
        path_control_cssplit (Union[Path, str]): Path of sample's jsonl containing CSSPLIT

    Returns:
        list[int]: _description_
    """
    # Create empty templates
    length = len(control_cssplit[0].split(","))
    sample_table = [[1, 1] for _ in range(length)]
    control_table = [[1, 1] for _ in range(length)]
    # Calculate match and mismatch
    for sample_cs, control_cs in zip(sample_cssplit, control_cssplit):
        sample_cs = sample_cs.split(",")
        control_cs = control_cs.split(",")
        for i, (s, c) in enumerate(zip(sample_cs, control_cs)):
            if s.startswith("="):
                sample_table[i][0] += 1
            else:
                sample_table[i][1] += 1
            if c.startswith("="):
                control_table[i][0] += 1
            else:
                control_table[i][1] += 1
    # Calculate match and mismatch
    different_loci = []
    for i, (s, c) in enumerate(zip(sample_table, control_table)):
        # odds = (s[1] * c[0]) / (s[0] * c[1])
        pval = st.chi2_contingency([s, c])[1]
        if pval < 0.01:  # and odds > 1
            different_loci.append(i)
    return different_loci


sample_table[3]
control_table[3]
