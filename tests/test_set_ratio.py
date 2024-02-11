from __future__ import annotations

import Levenshtein


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def test_documented():
    """
    run tests from documentation
    """
    ratio = Levenshtein.setratio(
        ["newspaper", "litter bin", "tinny", "antelope"],
        ["caribou", "sausage", "gorn", "woody"],
    )
    assert isclose(ratio, 0.2818452380952381)
