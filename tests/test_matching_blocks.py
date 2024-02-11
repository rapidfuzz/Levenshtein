from __future__ import annotations

from rapidfuzz.distance import MatchingBlock

from Levenshtein import editops, matching_blocks


def test_simple():
    a, b = "spam", "park"
    assert matching_blocks(editops(a, b), a, b) == [
        MatchingBlock(1, 0, 2),
        MatchingBlock(4, 4, 0),
    ]
    assert matching_blocks(editops(a, b), len(a), len(b)) == [
        MatchingBlock(1, 0, 2),
        MatchingBlock(4, 4, 0),
    ]
    assert matching_blocks(editops("", ""), 0, 0) == [MatchingBlock(0, 0, 0)]
    assert matching_blocks(editops("", "a"), 0, 1) == [MatchingBlock(0, 1, 0)]
    assert matching_blocks(editops("a", ""), 1, 0) == [MatchingBlock(1, 0, 0)]
    assert matching_blocks(editops("a", "a"), 1, 1) == [
        MatchingBlock(0, 0, 1),
        MatchingBlock(1, 1, 0),
    ]
