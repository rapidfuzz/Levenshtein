from typing import Callable, Hashable, List, Sequence, Tuple, Union, Optional, overload

__author__: str
__license__: str
__version__: str

_Editops = List[Tuple[str, int, int]]
_Opcodes = List[Tuple[str, int, int, int, int]]
_MatchingBlocks = List[Tuple[int, int, int]]
_AnyEditops = Union[_Editops, _Opcodes]

def inverse(edit_operations: list) -> list: ...

@overload
def editops(s1: str, s2: str) -> _Editops: ...

@overload
def editops(ops: _AnyEditops, s1: Union[str, int], s2: Union[str, int]) -> _Editops: ...

@overload
def opcodes(s1: str, s2: str) -> _Opcodes: ...

@overload
def opcodes(ops: _AnyEditops, s1: Union[str, int], s2: Union[str, int]) -> _Opcodes: ...

def matching_blocks(
    edit_operations: _AnyEditops,
    source_string: Union[str, int],
    destination_string: Union[str, int],
) -> _MatchingBlocks: ...
def subtract_edit(edit_operations: _Editops, subsequence: _Editops) -> _Editops: ...
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
    weights: Tuple[int, int, int] = (1, 1, 1),
    processor: Optional[Callable] = None,
    score_cutoff: Optional[float] = None,
) -> int: ...
def ratio(
    s1: Sequence[Hashable],
    s2: Sequence[Hashable],
    *,
    prefix_weight: float = 0.1,
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
    prefix_weight: float = 0.1,
    processor: Optional[Callable] = None,
    score_cutoff: Optional[float] = None,
) -> float: ...
