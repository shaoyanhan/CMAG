import pytest
from CMAG.pick_rep import calc_n50


def test_n50_success():
    assert 7 == calc_n50([10, 1, 9, 2, 8, 3, 7, 4, 6, 5])


def test_n50_failed():
    assert 10 == calc_n50([10, 1, 9, 2, 8, 3, 7, 4, 6, 5])