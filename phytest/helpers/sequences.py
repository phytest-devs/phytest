import re
from builtins import max as builtin_max
from typing import Optional
from warnings import warn

from Bio.SeqRecord import SeqRecord


def assert_sequence_valid_alphabet(sequence: SeqRecord, alphabet: str = "ATCGN-") -> None:
    regex_invalid = re.compile(f"[^{alphabet}]")
    assert not regex_invalid.search(str(sequence.seq)), f"Invalid pattern found in '{sequence.id}'."


def assert_sequence_length(
    sequence: SeqRecord,
    length: Optional[int] = None,
    *,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
) -> None:
    sequence_length = len(sequence.seq)
    if length is not None:
        assert sequence_length == length
    if min is not None:
        assert sequence_length >= min
    if max is not None:
        assert sequence_length <= max
    if warning is not None and sequence_length != warning:
        warn(f"Sequence length '{sequence_length}' != {warning}")


def assert_sequence_count(
    sequence: SeqRecord,
    pattern: str,
    *,
    count: Optional[int] = None,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
) -> None:
    base_count = sequence.seq.count(pattern)
    if count is not None:
        assert base_count == count
    if min is not None:
        assert base_count >= min
    if max is not None:
        assert base_count <= max
    if warning is not None and base_count != warning:
        warn(f"Count of '{pattern}' in {sequence.id} != {warning}")


def assert_sequence_count_Ns(
    sequence: SeqRecord,
    count: Optional[int] = None,
    *,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
) -> None:
    assert_sequence_count(sequence=sequence, pattern='N', count=count, min=min, max=max, warning=warning)


def assert_sequence_count_gaps(
    sequence: SeqRecord,
    count: Optional[int] = None,
    *,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
) -> None:
    assert_sequence_count(sequence=sequence, pattern='-', count=count, min=min, max=max, warning=warning)


def assert_sequence_longest_stretch(
    sequence: SeqRecord,
    pattern: str,
    *,
    count: Optional[int] = None,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
):
    matches = re.findall(f'{pattern}+', str(sequence.seq))
    longest_stretch = len(builtin_max(matches))
    if count is not None:
        assert longest_stretch == count
    if min is not None:
        assert longest_stretch >= min
    if max is not None:
        assert longest_stretch <= max
    if warning is not None and longest_stretch != warning:
        warn(f"Longest stretch of '{pattern}' in {sequence.id} != {warning}")


def assert_sequence_longest_stretch_Ns(
    sequence: SeqRecord,
    *,
    count: Optional[int] = None,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
):
    assert_sequence_longest_stretch(sequence=sequence, pattern='N', count=count, min=min, max=max, warning=warning)


def assert_sequence_longest_stretch_gaps(
    sequence: SeqRecord,
    *,
    count: Optional[int] = None,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
):
    assert_sequence_longest_stretch(sequence=sequence, pattern='-', count=count, min=min, max=max, warning=warning)


def assert_sequence_startswith(sequence: SeqRecord, pattern: str):
    assert sequence.seq.startswith(pattern)


def assert_sequence_endswith(sequence: SeqRecord, pattern: str):
    assert sequence.seq.endswith(pattern)


def assert_sequence_contains_motif(sequence: SeqRecord, motif: str):
    assert_sequence_count(sequence=sequence, pattern=motif, min=1)
