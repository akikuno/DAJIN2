# from pathlib import Path
# from src.DAJIN2.core.preprocess.get_index_mapping import calculate_index_mapping
# from src.DAJIN2.core.preprocess.get_index_mapping import get_index_mapping


# def test_calculate_index_mapping_substitution():
#     cssplits = ["=A", "=A", "*AG", "=A", "=A"]
#     test = calculate_index_mapping(cssplits)
#     answer = {0: 0, 1: 1, 3: 3, 4: 4}
#     assert test == answer


# def test_calculate_index_mapping_insertion():
#     cssplits = ["=A", "=A", "+G|+G|=A", "=A"]
#     test = calculate_index_mapping(cssplits)
#     answer = {0: 0, 1: 1, 2: 4, 3: 5}
#     assert test == answer


# def test_calculate_index_mapping_insertion_2():
#     cssplits = ["=A", "=A", "+G|+G|=A", "=A", "+G|+G|=A", "=A"]
#     test = calculate_index_mapping(cssplits)
#     answer = {0: 0, 1: 1, 2: 4, 3: 5, 4: 8, 5: 9}
#     assert test == answer


# def test_calculate_index_mapping_deletion():
#     cssplits = ["=A", "=A", "-G", "-G", "=A", "=A"]
#     test = calculate_index_mapping(cssplits)
#     answer = {0: 0, 1: 1, 4: 2, 5: 3}
#     assert test == answer


# def test_calculate_index_mapping_inversion():
#     cssplits = ["=A", "=g", "=g", "=A"]
#     test = calculate_index_mapping(cssplits)
#     answer = {0: 0, 1: 2, 2: 1, 3: 3}
#     assert test == answer


# def test_get_index_mapping_substitution():
#     test = get_index_mapping("tests/data/preprocess/get_index_mapping/substitution/")
#     assert len(test["albino"]) == 2844
#     assert 828 not in test["albino"]


# def test_get_index_mapping_insertion():
#     """
#     Cables2 flox design
#     - Control 1731 = Flox 1731
#     - Control 1732 = Flox 1778 due to left-loxp
#     - Control 2380 = Flox 2426 due to left-loxp
#     - Control 2381 = Flox 2489 due to both-loxp
#     """
#     test = get_index_mapping("tests/data/preprocess/get_index_mapping/insertion/")
#     assert test["flox"][1731] == 1731
#     assert test["flox"][1732] == 1778
#     assert test["flox"][2380] == 2426
#     assert test["flox"][2381] == 2489


# def test_get_index_mapping_deletion():
#     """
#     Stx2 2-cut deletion design
#     - deletion range: 2011~2737
#     """
#     test = get_index_mapping("tests/data/preprocess/get_index_mapping/deletion/")
#     assert 2010 in test["deletion"]
#     assert 2011 not in test["deletion"]
#     assert 2737 not in test["deletion"]
#     assert 2738 in test["deletion"]


# def test_get_index_mapping_inversion():
#     """
#     Stx2 2-cut deletion design
#     - inversion range: 2013~2739
#     """
#     test = get_index_mapping("tests/data/preprocess/get_index_mapping/inversion/")
#     assert test["inversion"][2012] == 2012
#     assert test["inversion"][2013] == 2739
#     assert test["inversion"][2739] == 2013
#     assert test["inversion"][2740] == 2740
