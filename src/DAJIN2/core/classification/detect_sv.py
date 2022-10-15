import re


def detect_sv(CSSPLIT: str, threshold: int = 50) -> bool:
    is_SV = False
    cssplit = CSSPLIT.replace(",", "")
    if "N" * threshold in cssplit:
        is_SV = True
    if re.search(rf"(\+[ACGTN]\|){{{threshold}}}", cssplit):
        is_SV = True
    if re.search(rf"(\-[ACGTN]){{{threshold}}}", cssplit):
        is_SV = True
    if re.search(rf"(\*[ACGTN][ACGTN]){{{threshold}}}", cssplit):
        is_SV = True
    if re.search(r"[acgtn]", cssplit):
        is_SV = True
    return is_SV
