import re
from builtins import max as builtin_max
from typing import List, Optional, Union

from Bio import AlignIO
from Bio import SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord

from ..utils import PhytestObject, assert_or_warn


class Sequence(PhytestObject, SeqRecord):
    @classmethod
    def parse(cls, alignment_path, alignment_format) -> 'Sequence':
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
            for r in SeqIO.parse(alignment_path, alignment_format)
        )

    def assert_valid_alphabet(self, alphabet: str = "ATCGN-", *, warning: bool = False) -> None:
        """
        Asserts that that the sequence only contains particular charaters.

        Args:
            alphabet (str): A string containing legal charaters. Defaults to 'ATCGN-'.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        regex_invalid = re.compile(f"[^{re.escape(alphabet)}]")
        result = regex_invalid.search(str(self.seq))
        if result:
            assert_or_warn(
                not result,
                warning,
                f"Invalid pattern found in '{self.id}'.",
                f"Character '{result.group(0)}' at position {result.start(0)+1} found which is not in alphabet '{alphabet}'.",
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
        Asserts that that the sequence length meets the specified criteria.

        Args:
            length (int, optional): If set, then sequence length must be equal to this value. Defaults to None.
            min (int, optional): If set, then sequence length must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then sequence length must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        sequence_length = len(self.seq)
        if length is not None:
            assert_or_warn(
                sequence_length == length,
                warning,
                f"Sequence length of '{self.id}' ({sequence_length}) is not equal to the required length of {length}.",
            )
        if min is not None:
            assert_or_warn(
                sequence_length >= min,
                warning,
                f"Sequence length of '{self.id}' ({sequence_length}) is less than the minimum {min}.",
            )
        if max is not None:
            assert_or_warn(
                sequence_length <= max,
                warning,
                f"Sequence length of '{self.id}' ({sequence_length}) is greater than the maximum {max}.",
            )

    def assert_count(
        self,
        pattern: str,
        *,
        count: Optional[int] = None,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ) -> None:
        """
        Asserts that the count of a pattern in the sequence meets the specified criteria.

        Args:
            pattern: (str): the pattern to count in the the sequence.
            count (int, optional): If set, then pattern count must be equal to this value. Defaults to None.
            min (int, optional): If set, then pattern count must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then pattern count must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        base_count = self.seq.count(pattern)
        summary = f"Sequence '{self.id}' matches pattern '{pattern}' {base_count} time(s)."
        if count is not None:
            assert_or_warn(
                base_count == count,
                warning,
                summary,
                f"This is not equal to the required number of {count}.",
            )
        if min is not None:
            assert_or_warn(
                base_count >= min,
                warning,
                summary,
                f"This is less than the minimum {min}.",
            )
        if max is not None:
            assert_or_warn(
                base_count <= max,
                warning,
                summary,
                f"This is greater than the maximum {max}.",
            )

    def assert_percent(
        self,
        nucleotide: Union[str, List[str]],
        *,
        percent: Optional[float] = None,
        min: Optional[float] = None,
        max: Optional[float] = None,
        warning: bool = False,
    ) -> None:
        """
        Asserts that the percentage of a nucleotide in the sequence meets the specified criteria.

        Args:
            nucleotide: (Union[str, List[str]]): The nucleotide(s) to count in the the sequence.
            percent (float, optional): If set, then nucleotide percentage must be equal to this value. Defaults to None.
            min (float, optional): If set, then nucleotide percentage must be equal to or greater than this value. Defaults to None.
            max (float, optional): If set, then nucleotide percentage must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        try:
            if isinstance(nucleotide, str):
                if len(nucleotide) > 1:
                    raise ValueError(
                        f"The length of the requested nucleotide '{nucleotide}' is more than a single character. "
                        "This value should either be a single character (i.e. A, G, C, T) or a list of single characters."
                    )
                base_percent = (self.seq.count(nucleotide) * 100.0) / len(self.seq)
            elif isinstance(nucleotide, list):
                base_percent = (sum(self.seq.count(x) for x in nucleotide) * 100) / len(self.seq)
                nucleotide = ', '.join(nucleotide)
            else:
                raise ValueError(f"Nucleotide must be str or list and cannot be of type '{type(nucleotide)}'.")
        except ZeroDivisionError:
            base_percent = 0.0
        summary = f"Sequence '{self.id}' contains {base_percent} percent '{nucleotide}'."
        if percent is not None:
            assert_or_warn(
                base_percent == percent,
                warning,
                summary,
                f"This is not equal to the required percentage of {percent}.",
            )
        if min is not None:
            assert_or_warn(
                base_percent >= min,
                warning,
                summary,
                f"This is less than the minimum {min}.",
            )
        if max is not None:
            assert_or_warn(
                base_percent <= max,
                warning,
                summary,
                f"This is greater than the maximum {max}.",
            )

    def assert_percent_GC(
        self,
        percent: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ) -> None:
        """
        Asserts that the percent of GC's (ambiguous nucleotide S) in the sequence meets the specified criteria.

        Args:
            percent (float, optional): If set, then the percentage of GC's must be equal to this value. Defaults to None.
            min (float, optional): If set, then the percentage of GC's must be equal to or greater than this value. Defaults to None.
            max (float, optional): If set, then the percentage of GC's must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        self.assert_percent(
            nucleotide=["G", "C", "g", "c", "S", "s"], percent=percent, min=min, max=max, warning=warning
        )

    def assert_percent_N(
        self,
        percent: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ) -> None:
        """
        Asserts that the percent of N's in the sequence meets the specified criteria.

        Args:
            percent (float, optional): If set, then the percentage of N's must be equal to this value. Defaults to None.
            min (float, optional): If set, then the percentage of N's must be equal to or greater than this value. Defaults to None.
            max (float, optional): If set, then the percentage of N's must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        self.assert_percent(nucleotide=["N", "n"], percent=percent, min=min, max=max, warning=warning)

    def assert_percent_gaps(
        self,
        percent: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ) -> None:
        """
        Asserts that the percent of gaps (-) in the sequence meets the specified criteria.

        Args:
            percent (float, optional): If set, then the percentage of gaps must be equal to this value. Defaults to None.
            min (float, optional): If set, then the percentage of gaps must be equal to or greater than this value. Defaults to None.
            max (float, optional): If set, then the percentage of gaps must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        self.assert_percent(nucleotide='-', percent=percent, min=min, max=max, warning=warning)

    def assert_count_Ns(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ) -> None:
        """
        Asserts that the number of a N's in the sequence meets the specified criteria.

        Args:
            count (int, optional): If set, then the number of N's must be equal to this value. Defaults to None.
            min (int, optional): If set, then the number of N's must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then the number of N's must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        self.assert_count(pattern='N', count=count, min=min, max=max, warning=warning)

    def assert_count_gaps(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ) -> None:
        """
        Asserts that the number of a gaps (-) in the sequence meets the specified criteria.

        Args:
            count (int, optional): If set, then the number of gaps (-) must be equal to this value. Defaults to None.
            min (int, optional): If set, then the number of gaps (-) must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then the number of gaps (-) must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        self.assert_count(pattern='-', count=count, min=min, max=max, warning=warning)

    def assert_longest_stretch(
        self,
        pattern: str,
        *,
        count: Optional[int] = None,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ):
        """
        Asserts that the longest stretch of a pattern in the sequence meets the specified criteria.

        e.g. the longest stretch of N's in 'ANNNANNA' is 3.

        Args:
            pattern: (str): the pattern to count in the the sequence.
            count (int, optional): If set, then the longest stretch of the pattern must be equal to this value. Defaults to None.
            min (int, optional): If set, then the longest stretch of the pattern must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then the longest stretch of the pattern must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        matches = re.findall(f'{pattern}+', str(self.seq))
        longest_stretch = len(builtin_max(matches)) if matches else 0
        summary = f"The longest stretch of pattern '{pattern}' in sequence '{self.id}' is {longest_stretch}."
        if count is not None:
            assert_or_warn(
                longest_stretch == count,
                warning,
                summary,
                f"This is not equal to the required number of {count}.",
            )
        if min is not None:
            assert_or_warn(
                longest_stretch >= min,
                warning,
                summary,
                f"This is less than the minimum {min}.",
            )
        if max is not None:
            assert_or_warn(
                longest_stretch <= max,
                warning,
                summary,
                f"This is greater than the maximum {max}.",
            )

    def assert_longest_stretch_Ns(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ):
        """
        Asserts that the longest stretch of a N's in the sequence meets the specified criteria.

        e.g. the logest stretch of N's in 'ANNNANNA' is 3.

        Args:
            count (int, optional): If set, then the longest stretch of N's must be equal to this value. Defaults to None.
            min (int, optional): If set, then the longest stretch of N's must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then the longest stretch of N's must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        self.assert_longest_stretch(pattern='N', count=count, min=min, max=max, warning=warning)

    def assert_longest_stretch_gaps(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: bool = False,
    ):
        """
        Asserts that the longest stretch of a gaps (-) in the sequence meets the specified criteria.

        e.g. the logest stretch of gaps (-) in 'A---A--A' is 3.

        Args:
            count (int, optional): If set, then the longest stretch of gaps (-) must be equal to this value. Defaults to None.
            min (int, optional): If set, then the longest stretch of gaps (-) must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then the longest stretch of gaps (-) must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        self.assert_longest_stretch(pattern='-', count=count, min=min, max=max, warning=warning)

    def assert_startswith(self, pattern: str, *, warning: bool = False):
        """
        Asserts that the sequence starts with a particular pattern.

        Args:
            pattern (str): The sequence must start with this value.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        assert_or_warn(
            self.seq.startswith(pattern),
            warning,
            f"Sequence '{self.id}' does not start with '{pattern}'.",
        )

    def assert_endswith(self, pattern: str, *, warning: bool = False):
        """
        Asserts that the sequence ends with a particular pattern.

        Args:
            pattern (str): The sequence must end with this value.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        assert_or_warn(
            self.seq.endswith(pattern),
            warning,
            f"Sequence '{self.id}' does not end with '{pattern}'.",
        )

    def assert_contains(self, pattern: str, *, warning: bool = False):
        """
        Asserts that the sequence contains a particular pattern.

        Args:
            pattern (str): The sequence must contain this value.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        self.assert_count(pattern=pattern, min=1, warning=warning)
