from src.DAJIN2.consensus import module_consensus
from importlib import reload

reload(module_consensus)


def test_module_consensus():
    cssplits = ["=A,*GA,=C", "=A,*GA,-C", "=A,*GA,-C"]
    test = module_consensus.call(cssplits)
    answer = "=A,*GA,-C"
    assert test == answer

