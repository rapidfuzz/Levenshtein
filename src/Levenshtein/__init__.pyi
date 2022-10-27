from typing import Callable, Hashable, List, Sequence, Tuple, Union, Optional, overload

__author__: str
__license__: str
__version__: str

_EditopsList = List[Tuple[str, int, int]]
_OpcodesList = List[Tuple[str, int, int, int, int]]
_MatchingBlocks = List[Tuple[int, int, int]]
_AnyEditops = Union[_EditopsList, _OpcodesList]

def inverse(edit_operations: list) -> list: ...
@overload
def editops(s1: Sequence[Hashable], s2: Sequence[Hashable]) -> _EditopsList: ...
@overload
def editops(
    ops: _AnyEditops,
    s1: Union[Sequence[Hashable], int],
    s2: Union[Sequence[Hashable], int],
) -> _EditopsList: ...
@overload
def opcodes(s1: Sequence[Hashable], s2: Sequence[Hashable]) -> _OpcodesList: ...
@overload
def opcodes(
    ops: _AnyEditops,
    s1: Union[Sequence[Hashable], int],
    s2: Union[Sequence[Hashable], int],
) -> _OpcodesList: ...
def matching_blocks(
    edit_operations: _AnyEditops,
    source_string: Union[Sequence[Hashable], int],
    destination_string: Union[Sequence[Hashable], int],
) -> _MatchingBlocks: ...
def subtract_edit(
    edit_operations: _EditopsList, subsequence: _EditopsList
) -> _EditopsList: ...
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
