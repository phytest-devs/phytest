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
        if format in AlignIO._FormatToIterator:
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
        else:
            raise ValueError("Unknown format '%s'" % format)

    def assert_valid_alphabet(self, alphabet: str = "ATCGN-") -> None:
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
        self.assert_count(pattern='N', count=count, min=min, max=max, warning=warning)

    def assert_count_gaps(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ) -> None:
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
        self.assert_longest_stretch(pattern='N', count=count, min=min, max=max, warning=warning)

    def assert_longest_stretch_gaps(
        self,
        count: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ):
        self.assert_longest_stretch(pattern='-', count=count, min=min, max=max, warning=warning)

    def assert_startswith(self, pattern: str):
        assert self.seq.startswith(pattern)

    def assert_endswith(self, pattern: str):
        assert self.seq.endswith(pattern)

    def assert_contains_motif(self, motif: str):
        self.assert_count(pattern=motif, min=1)
