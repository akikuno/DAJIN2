from DAJIN2.core.preprocess.midsv_caller import replace_n_to_d


def test_replace_n_to_d():
    midsv_sample = [{"CSSPLIT": "N,N,N,=A,N,=C,N,N"}, {"CSSPLIT": "N,N,=A,N,N,=C,=C,N"}]
    sequence = "GCAACCCC"
    test = replace_n_to_d(midsv_sample, sequence)
    test = list(test)
    answer = [{"CSSPLIT": "N,N,N,=A,-C,=C,N,N"}, {"CSSPLIT": "N,N,=A,-A,-C,=C,=C,N"}]
    assert test == answer


def test_replace_n_to_d_large_N():
    midsv_sample = [{"CSSPLIT": "N,N,N,N,N,N,=C,N,=A"}]
    sequence = "GCAACCCCA"
    test = replace_n_to_d(midsv_sample, sequence)
    test = list(test)
    answer = [{"CSSPLIT": "N,N,N,N,N,N,=C,-C,=A"}]
    assert test == answer
