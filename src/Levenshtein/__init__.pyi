from typing import Callable, Hashable, List, Sequence, Tuple, Union, Optional, overload
from rapidfuzz.distance import Editops, Opcodes, MatchingBlocks

__author__: str
__license__: str
__version__: str

_AnyEditops = Union[Editops, Opcodes]

def inverse(edit_operations: list) -> list: ...
@overload
def editops(s1: str, s2: str) -> Editops: ...
@overload
def editops(
    ops: _AnyEditops, s1: Union[str, int], s2: Union[str, int]
) -> Editops: ...
@overload
def opcodes(s1: str, s2: str) -> Opcodes: ...
@overload
def opcodes(
    ops: _AnyEditops, s1: Union[str, int], s2: Union[str, int]
) -> Opcodes: ...
def matching_blocks(
    edit_operations: _AnyEditops,
    source_string: Union[str, int],
    destination_string: Union[str, int],
) -> MatchingBlocks: ...
def subtract_edit(
    edit_operations: Editops, subsequence: Editops
) -> Editops: ...
def apply_edit(
    edit_operations: _AnyEditops, source_string: str, destination_string: str
) -> str: ...
def median(
    strlist: List[Union[str, bytes]], wlist: Optional[List[float]] = None
) -> str: ...
def quickmedian(
    strlist: List[Union[str, bytes]], wlist: Optional[List[float]] = None
) -> str: ...
def median_improve(
    string: Union[str, bytes],
    strlist: List[Union[str, bytes]],
    wlist: Optional[List[float]] = None,
) -> str: ...
def setmedian(
    strlist: List[Union[str, bytes]], wlist: Optional[List[float]] = None
) -> str: ...
def setratio(
    strlist1: List[Union[str, bytes]], strlist2: List[Union[str, bytes]]
) -> float: ...
def seqratio(
    strlist1: List[Union[str, bytes]], strlist2: List[Union[str, bytes]]
) -> float: ...
def distance(
    s1: Sequence[Hashable],
    s2: Sequence[Hashable],
    *,
    weights: Optional[Tuple[int, int, int]] = (1, 1, 1),
    processor: Optional[Callable] = None,
    score_cutoff: Optional[float] = None,
) -> int: ...
def ratio(
    s1: Sequence[Hashable],
    s2: Sequence[Hashable],
    *,
    processor: Optional[Callable] = None,
    score_cutoff: Optional[float] = None,
) -> float: ...
def hamming(
    s1: Sequence[Hashable],
    s2: Sequence[Hashable],
    *,
    processor: Optional[Callable] = None,
    score_cutoff: Optional[float] = None,
) -> int: ...
def jaro(
    s1: Sequence[Hashable],
    s2: Sequence[Hashable],
    *,
    processor: Optional[Callable] = None,
    score_cutoff: Optional[float] = None,
) -> float: ...
def jaro_winkler(
    s1: Sequence[Hashable],
    s2: Sequence[Hashable],
    *,
    prefix_weight: Optional[float] = 0.1,
    processor: Optional[Callable] = None,
    score_cutoff: Optional[float] = None,
) -> float: ...
