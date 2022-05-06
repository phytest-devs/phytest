import re
from builtins import max as builtin_max
from typing import Optional
from warnings import warn

from Bio import AlignIO
from Bio import SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord


class Sequence(SeqRecord):
    @classmethod
    def parse(cls, alignment_path, alignment_format) -> 'Sequence':
        # Use Bio.AlignIO to read in the alignments
        return (
            Sequence(
                r.seq,
                id=r.id,
                name=r.name,
                description=r.description,
                dbxrefs=r.dbxrefs,
                features=r.features,
                annotations=r.annotations,
                letter_annotations=r.letter_annotations,
            )
            for alignment in AlignIO.parse(alignment_path, alignment_format)
            for r in alignment
        )

    def assert_valid_alphabet(self, alphabet: str = "ATCGN-") -> None:
        """
        Asserts that that the sequence only contains particular charaters.

        Args:
            alphabet (str): A string containing legal charaters. Defaults to 'ATCGN-'.
        """
        regex_invalid = re.compile(f"[^{alphabet}]")
        assert not regex_invalid.search(str(self.seq)), f"Invalid pattern found in '{self.id}'!"

    def assert_length(
        self,
        length: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ) -> None:
        """
        Asserts that that the sequence length meets the specified criteria.

        Args:
            length (Optional[int], optional): If set, then sequence length must be equal to this value. Defaults to None.
            min (Optional[int], optional): If set, then sequence length must be equal to or greater than this value. Defaults to None.
            max (Optional[int], optional): If set, then sequence length must be equal to or less than this value. Defaults to None.
            warning (Optional[int], optional): If set, raise a warning if the sequence length is not equal to this value. Defaults to None.
        """
        sequence_length = len(self.seq)
        if length is not None:
            assert sequence_length == length
        if min is not None:
            assert sequence_length >= min
        if max is not None:
            assert sequence_length <= max
        if warning is not None and sequence_length != warning:
            warn(f"Sequence length '{sequence_length}' != {warning}")

    def assert_count(
        self,
        pattern: str,
        *,
        count: Optional[int] = None,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ) -> None:
        """
        Asserts that the count of a pattern in the sequence meets the specified criteria.

        Args:
            pattern: (str): the pattern to count in the the sequence.
            count (Optional[int], optional): If set, then pattern count must be equal to this value. Defaults to None.
            min (Optional[int], optional): If set, then pattern count must be equal to or greater than this value. Defaults to None.
            max (Optional[int], optional): If set, then pattern count must be equal to or less than this value. Defaults to None.
            warning (Optional[int], optional): If set, raise a warning if the pattern count is not equal to this value. Defaults to None.
        """
        base_count = self.seq.count(pattern)
        if count is not None:
            assert base_count == count
        if min is not None:
            assert base_count >= min
        if max is not None:
            assert base_count <= max
        if warning is not None and base_count != warning:
            warn(f"Count of '{pattern}' in {self.id} != {warning}")

    def assert_count_Ns(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ) -> None:
        """
        Asserts that the number of a N's in the sequence meets the specified criteria.

        Args:
            count (Optional[int], optional): If set, then the number of N's must be equal to this value. Defaults to None.
            min (Optional[int], optional): If set, then the number of N's must be equal to or greater than this value. Defaults to None.
            max (Optional[int], optional): If set, then the number of N's must be equal to or less than this value. Defaults to None.
            warning (Optional[int], optional): If set, raise a warning if the number of N's is not equal to this value. Defaults to None.
        """
        self.assert_count(pattern='N', count=count, min=min, max=max, warning=warning)

    def assert_count_gaps(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ) -> None:
        """
        Asserts that the number of a gaps (-) in the sequence meets the specified criteria.

        Args:
            count (Optional[int], optional): If set, then the number of gaps (-) must be equal to this value. Defaults to None.
            min (Optional[int], optional): If set, then the number of gaps (-) must be equal to or greater than this value. Defaults to None.
            max (Optional[int], optional): If set, then the number of gaps (-) must be equal to or less than this value. Defaults to None.
            warning (Optional[int], optional): If set, raise a warning if the number of gaps (-) is not equal to this value. Defaults to None.
        """
        self.assert_count(pattern='-', count=count, min=min, max=max, warning=warning)

    def assert_longest_stretch(
        self,
        pattern: str,
        *,
        count: Optional[int] = None,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ):
        """
        Asserts that the longest stretch of a pattern in the sequence meets the specified criteria.

        e.g. the logest stretch of N's in 'ANNNANNA' is 3.

        Args:
            pattern: (str): the pattern to count in the the sequence.
            count (Optional[int], optional): If set, then the longest stretch of the pattern must be equal to this value. Defaults to None.
            min (Optional[int], optional): If set, then the longest stretch of the pattern must be equal to or greater than this value. Defaults to None.
            max (Optional[int], optional): If set, then the longest stretch of the pattern must be equal to or less than this value. Defaults to None.
            warning (Optional[int], optional): If set, raise a warning if the longest stretch of the pattern is not equal to this value. Defaults to None.
        """
        matches = re.findall(f'{pattern}+', str(self.seq))
        if matches:
            longest_stretch = len(builtin_max(matches))
        else:
            longest_stretch = 0
        if count is not None:
            assert longest_stretch == count
        if min is not None:
            assert longest_stretch >= min
        if max is not None:
            assert longest_stretch <= max, f"Longest stretch of '{pattern}' in '{self.id}' > {max}!"
        if warning is not None and longest_stretch != warning:
            warn(f"Longest stretch of '{pattern}' in {self.id} != {warning}")

    def assert_longest_stretch_Ns(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ):
        """
        Asserts that the longest stretch of a N's in the sequence meets the specified criteria.

        e.g. the logest stretch of N's in 'ANNNANNA' is 3.

        Args:
            count (Optional[int], optional): If set, then the longest stretch of N's must be equal to this value. Defaults to None.
            min (Optional[int], optional): If set, then the longest stretch of N's must be equal to or greater than this value. Defaults to None.
            max (Optional[int], optional): If set, then the longest stretch of N's must be equal to or less than this value. Defaults to None.
            warning (Optional[int], optional): If set, raise a warning if the longest stretch of N's is not equal to this value. Defaults to None.
        """
        self.assert_longest_stretch(pattern='N', count=count, min=min, max=max, warning=warning)

    def assert_longest_stretch_gaps(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ):
        """
        Asserts that the longest stretch of a gaps (-) in the sequence meets the specified criteria.

        e.g. the logest stretch of gaps (-) in 'A---A--A' is 3.

        Args:
            count (Optional[int], optional): If set, then the longest stretch of gaps (-) must be equal to this value. Defaults to None.
            min (Optional[int], optional): If set, then the longest stretch of gaps (-) must be equal to or greater than this value. Defaults to None.
            max (Optional[int], optional): If set, then the longest stretch of gaps (-) must be equal to or less than this value. Defaults to None.
            warning (Optional[int], optional): If set, raise a warning if the longest stretch of gaps (-) is not equal to this value. Defaults to None.
        """
        self.assert_longest_stretch(pattern='-', count=count, min=min, max=max, warning=warning)

    def assert_startswith(self, pattern: str):
        """
        Asserts that the sequence starts with a particular pattern.

        Args:
            pattern (str): The sequence must start with this value.
        """
        assert self.seq.startswith(pattern)

    def assert_endswith(self, pattern: str):
        """
        Asserts that the sequence ends with a particular pattern.

        Args:
            pattern (str): The sequence must end with this value.
        """
        assert self.seq.endswith(pattern)

    def assert_contains(self, pattern: str):
        """
        Asserts that the sequence contains a particular pattern.

        Args:
            pattern (str): The sequence must contain this value.
        """
        self.assert_count(pattern=pattern, min=1)
