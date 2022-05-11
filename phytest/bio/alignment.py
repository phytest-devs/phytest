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
        """
        Asserts that the alignment width (the number of bases in the sequences) meets the specified criteria.

        Args:
            length (Optional[int], optional): If set, then alignment width must be equal to this value. Defaults to None.
            min (Optional[int], optional): If set, then alignment width must be equal to or greater than this value. Defaults to None.
            max (Optional[int], optional): If set, then alignment width must be equal to or less than this value. Defaults to None.
            warning (Optional[int], optional): If set, raise a warning if the alignment width is not equal to this value. Defaults to None.
        """
        alignment_width = self.get_alignment_length()
        if width is not None:
            assert alignment_width == width, f"Alignment width '{alignment_width}' != '{width}'"
        if min is not None:
            assert alignment_width >= min, f"Alignment width '{alignment_width}' < '{min}'"
        if max is not None:
            assert alignment_width <= max, f"Alignment width '{alignment_width}' > '{max}'"
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
        """
        Asserts that the alignment length (the number of sequences in the alignment) meets the specified criteria.

        Args:
            length (Optional[int], optional): If set, then alignment length must be equal to this value. Defaults to None.
            min (Optional[int], optional): If set, then alignment length must be equal to or greater than this value. Defaults to None.
            max (Optional[int], optional): If set, then alignment length must be equal to or less than this value. Defaults to None.
            warning (Optional[int], optional): If set, raise a warning if the alignment length is not equal to this value. Defaults to None.
        """
        alignment_length = len(self)
        if length is not None:
            assert alignment_length == length, f"Alignment length '{alignment_length}' != '{length}'"
        if min is not None:
            assert alignment_length >= min, f"Alignment length '{alignment_length}' < '{min}'"
        if max is not None:
            assert alignment_length <= max, f"Alignment length '{alignment_length}' > '{max}'"
        if warning is not None and alignment_length != warning:
            warn(f"Alignment length '{alignment_length}' != {warning}")
