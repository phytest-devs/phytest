from typing import Optional
from warnings import warn

from Bio.Align import MultipleSeqAlignment


def assert_alignment_width(
    alignment: MultipleSeqAlignment,
    width: Optional[int] = None,
    *,
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


def assert_alignment_length(
    alignment: MultipleSeqAlignment,
    length: Optional[int] = None,
    *,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
) -> None:
    alignment_length = len(alignment)
    if length is not None:
        assert alignment_length == length
    if min is not None:
        assert alignment_length > min
    if max is not None:
        assert alignment_length < max
    if warning is not None and alignment_length != warning:
        warn(f"Alignment length '{alignment_length}' != {warning}")
