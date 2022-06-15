from typing import Optional
from warnings import warn

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from ..utils import PhytestObject, assert_or_warn


class Alignment(PhytestObject, MultipleSeqAlignment):
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
        warning: bool = False,
    ) -> None:
        """
        Asserts that the alignment width (the number of bases in the sequences) meets the specified criteria.

        Args:
            length (int, optional): If set, then alignment width must be equal to this value. Defaults to None.
            min (int, optional): If set, then alignment width must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then alignment width must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        alignment_width = self.get_alignment_length()
        summary = f"The width of the alignment is {alignment_width}."

        if width is not None:
            assert_or_warn(
                alignment_width == width,
                warning,
                summary,
                f"This is not equal to the required width of {width}.",
            )
        if min is not None:
            assert_or_warn(
                alignment_width >= min,
                warning,
                summary,
                f"This is less than the minimum width of {min}.",
            )
        if max is not None:
            assert_or_warn(
                alignment_width <= max,
                warning,
                summary,
                f"This is greater than the maximum width of {max}.",
            )

    def assert_length(
        self,
        length: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ) -> None:
        """
        Asserts that the alignment length (the number of sequences in the alignment) meets the specified criteria.

        Args:
            length (int, optional): If set, then alignment length must be equal to this value. Defaults to None.
            min (int, optional): If set, then alignment length must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then alignment length must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        alignment_length = len(self)
        summary = f"The number of sequences in the alignment is {alignment_length}."

        if length is not None:
            assert_or_warn(
                alignment_length == length,
                warning,
                summary,
                f"This is less than required number of {length}.",
            )
        if min is not None:
            assert_or_warn(
                alignment_length >= min,
                warning,
                summary,
                f"This is less than the minimum {min}.",
            )
        if max is not None:
            assert_or_warn(
                alignment_length <= max,
                warning,
                summary,
                f"This is greater than the maximum {max}.",
            )
