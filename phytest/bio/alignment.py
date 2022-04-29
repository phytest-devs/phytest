from typing import Optional
from warnings import warn

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


class Alignment(MultipleSeqAlignment):
    @classmethod
    def read(cls, alignment_path, alignment_format) -> 'Alignment':
        alignment = AlignIO.read(alignment_path, alignment_format)
        return Alignment(
            alignment._records, annotations=alignment.annotations, column_annotations=alignment.column_annotations
        )

    def assert_width(
        self,
        width: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ) -> None:
        alignment_width = self.get_alignment_length()
        if width is not None:
            assert alignment_width == width
        if min is not None:
            assert alignment_width >= min
        if max is not None:
            assert alignment_width <= max
        if warning is not None and alignment_width != warning:
            warn(f"Alignment width '{alignment_width}' != {warning}")

    def assert_length(
        self,
        length: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ) -> None:
        alignment_length = len(self)
        if length is not None:
            assert alignment_length == length
        if min is not None:
            assert alignment_length >= min
        if max is not None:
            assert alignment_length <= max
        if warning is not None and alignment_length != warning:
            warn(f"Alignment length '{alignment_length}' != {warning}")
