import re
from typing import Optional
from warnings import warn

from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord


def assert_valid_alphabet(sequence: SeqRecord, *, alphabet: str = "ATCGN-") -> None:
    regex_invalid = re.compile(f"[^{alphabet}]")
    assert not regex_invalid.search(sequence.seq), f"Invalid base found in '{sequence.id}'."


def assert_length(
    sequence: SeqRecord,
    *,
    length: Optional[int] = None,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
) -> None:
    sequence_length = len(sequence.seq)
    if length is not None:
        assert sequence_length == length
    if min is not None:
        assert sequence_length > min
    if max is not None:
        assert sequence_length < max
    if warning is not None and sequence_length != warning:
        warn(f"Sequence length '{sequence_length}' != {warning}")


def assert_count(
    sequence: SeqRecord,
    *,
    base: str,
    count: Optional[int] = None,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
) -> None:
    base_count = sequence.seq.count(base)
    if count is not None:
        assert base_count == count
    if min is not None:
        assert count > min
    if max is not None:
        assert count < max
    if warning is not None and count != warning:
        warn(f"Count of '{base}' in {sequence.id} > {warning}")


def assert_alignment_width(
    alignment: MultipleSeqAlignment,
    *,
    width: Optional[int] = None,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
) -> None:
    alignment_width = alignment.get_alignment_length()
    if width is not None:
        assert alignment_width == width
    if min is not None:
        assert alignment_width > min
    if max is not None:
        assert alignment_width < max
    if warning is not None and alignment_width != warning:
        warn(f"Alignment width '{alignment_width}' != {warning}")


def assert_number_of_records(
    alignment: MultipleSeqAlignment,
    *,
    records: Optional[int] = None,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
) -> None:
    alignment_length = len(alignment)
    if records is not None:
        assert alignment_length == records
    if min is not None:
        assert alignment_length > min
    if max is not None:
        assert alignment_length < max
    if warning is not None and alignment_length != warning:
        warn(f"Alignment length '{alignment_length}' != {warning}")
