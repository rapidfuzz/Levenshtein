from __future__ import annotations

import Levenshtein


def test_weight_zero():
    """
    Strings with zero weights should be ignored
    """
    assert Levenshtein.quickmedian(["tes", "teste"], [0, 0]) == ""
    assert Levenshtein.quickmedian(["tes", "teste"], [1, 0]) == "tes"
    assert Levenshtein.quickmedian(["tes", "teste"], [0, 1]) == "teste"


def test_documented():
    """
    run tests from documentation
    """
    assert Levenshtein.median(["SpSm", "mpamm", "Spam", "Spa", "Sua", "hSam"]) == "Spam"
    fixme = [
        "Levnhtein",
        "Leveshein",
        "Leenshten",
        "Leveshtei",
        "Lenshtein",
        "Lvenstein",
        "Levenhtin",
        "evenshtei",
    ]
    assert Levenshtein.median(fixme) == "Levenshtein"
    assert Levenshtein.quickmedian(fixme) == "Levnshein"
    assert Levenshtein.median_improve("spam", fixme) == "enhtein"
    assert Levenshtein.median_improve(Levenshtein.median_improve("spam", fixme), fixme) == "Levenshtein"
    assert (
        Levenshtein.setmedian(
            [
                "ehee",
                "cceaes",
                "chees",
                "chreesc",
                "chees",
                "cheesee",
                "cseese",
                "chetese",
            ]
        )
        == "chees"
    )
