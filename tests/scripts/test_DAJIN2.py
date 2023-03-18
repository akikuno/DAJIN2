from src.DAJIN2 import DAJIN2
import os

def test_threads_ok():
    threads = 2
    test = DAJIN2.update_threads(threads)
    assert test == 2

def test_threads_minus():
    threads = -100
    test = DAJIN2.update_threads(threads)
    assert test == 1

def test_threads_over():
    threads = 10**100
    test = DAJIN2.update_threads(threads)
    assert test == os.cpu_count() - 1
