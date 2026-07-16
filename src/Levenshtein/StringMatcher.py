from __future__ import annotations

from typing import TYPE_CHECKING, Literal
from warnings import warn

from Levenshtein import distance, editops, matching_blocks, opcodes, ratio

if TYPE_CHECKING:
    from Levenshtein import _EditopsList, _MatchingBlocks, _OpcodesList


class StringMatcher:
    """A SequenceMatcher-like class built on the top of Levenshtein"""

    _ratio: float | None
    _distance: int | None
    _opcodes: _OpcodesList | None
    _editops: _EditopsList | None
    _matching_blocks: _MatchingBlocks | None

    def _reset_cache(self) -> None:
        self._ratio = self._distance = None
        self._opcodes = self._editops = self._matching_blocks = None

    def __init__(
        self,
        isjunk: Literal[False] | None = None,
        seq1: str = "",
        seq2: str = "",
        autojunk: Literal[False] = False,
    ) -> None:
        if isjunk:
            warn("isjunk NOT implemented, it will be ignored", stacklevel=1)
        if autojunk:
            warn("autojunk NOT implemented, it will be ignored", stacklevel=1)
        self._str1, self._str2 = seq1, seq2
        self._reset_cache()

    def set_seqs(self, seq1: str, seq2: str) -> None:
        self._str1, self._str2 = seq1, seq2
        self._reset_cache()

    def set_seq1(self, seq1: str) -> None:
        self._str1 = seq1
        self._reset_cache()

    def set_seq2(self, seq2: str) -> None:
        self._str2 = seq2
        self._reset_cache()

    def get_opcodes(self) -> _OpcodesList:
        if not self._opcodes:
            if self._editops:
                self._opcodes = opcodes(self._editops, self._str1, self._str2)
            else:
                self._opcodes = opcodes(self._str1, self._str2)
        return self._opcodes

    def get_editops(self) -> _EditopsList:
        if not self._editops:
            if self._opcodes:
                self._editops = editops(self._opcodes, self._str1, self._str2)
            else:
                self._editops = editops(self._str1, self._str2)
        return self._editops

    def get_matching_blocks(self) -> _MatchingBlocks:
        if not self._matching_blocks:
            self._matching_blocks = matching_blocks(self.get_opcodes(), self._str1, self._str2)
        return self._matching_blocks

    def ratio(self) -> float:
        if not self._ratio:
            self._ratio = ratio(self._str1, self._str2)
        return self._ratio

    def quick_ratio(self) -> float:
        # This is usually quick enough :o)
        if not self._ratio:
            self._ratio = ratio(self._str1, self._str2)
        return self._ratio

    def real_quick_ratio(self) -> float:
        len1, len2 = len(self._str1), len(self._str2)
        return 2.0 * min(len1, len2) / (len1 + len2)

    def distance(self) -> int:
        if not self._distance:
            self._distance = distance(self._str1, self._str2)
        return self._distance
