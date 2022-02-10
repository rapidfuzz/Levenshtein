import Levenshtein

def test_weight_zero():
    """
    Strings with zero weights should be ignored
    """
    assert Levenshtein.quickmedian(["tes", "teste"], [0, 0]) == ""
    assert Levenshtein.quickmedian(["tes", "teste"], [1,0]) == "tes"
    assert Levenshtein.quickmedian(["tes", "teste"], [0,1]) == "teste"
