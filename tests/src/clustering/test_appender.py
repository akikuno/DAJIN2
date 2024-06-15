from src.DAJIN2.core.clustering.appender import add_percent, add_readnum


def test_add_readnum():
    sample_data = [
        {"LABEL": 1},
        {"LABEL": 1},
        {"LABEL": 2},
        {"LABEL": 3},
        {"LABEL": 3},
        {"LABEL": 3},
    ]
    expected_output = [
        {"LABEL": 1, "READNUM": 2},
        {"LABEL": 1, "READNUM": 2},
        {"LABEL": 2, "READNUM": 1},
        {"LABEL": 3, "READNUM": 3},
        {"LABEL": 3, "READNUM": 3},
        {"LABEL": 3, "READNUM": 3},
    ]
    assert add_readnum(sample_data) == expected_output


def test_add_percent():
    sample_data = [
        {"LABEL": 1},
        {"LABEL": 1},
        {"LABEL": 2},
        {"LABEL": 3},
        {"LABEL": 3},
        {"LABEL": 3},
    ]
    # 1 occurs 2/6 times, 2 occurs 1/6 times, 3 occurs 3/6 times
    # percentages: 1 -> 33.333, 2 -> 16.667, 3 -> 50.000
    expected_output = [
        {"LABEL": 1, "PERCENT": 33.333},
        {"LABEL": 1, "PERCENT": 33.333},
        {"LABEL": 2, "PERCENT": 16.667},
        {"LABEL": 3, "PERCENT": 50.0},
        {"LABEL": 3, "PERCENT": 50.0},
        {"LABEL": 3, "PERCENT": 50.0},
    ]
    assert add_percent(sample_data) == expected_output
