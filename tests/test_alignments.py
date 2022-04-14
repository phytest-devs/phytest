import pytest
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phytest import alignments


def test_assert_alignment_width():
    alignment_path = 'examples/data/invalid.fasta'
    alignment = AlignIO.read(alignment_path, 'fasta')
    alignments.assert_alignment_width(alignment, width=100)
    alignments.assert_alignment_width(alignment, min=99)
    alignments.assert_alignment_width(alignment, max=101)
    with pytest.raises(AssertionError):
        alignments.assert_alignment_width(alignment, width=99)
        alignments.assert_alignment_width(alignment, min=101)
        alignments.assert_alignment_width(alignment, max=99)


def test_assert_alignment_length():
    alignment_path = 'examples/data/invalid.fasta'
    alignment = AlignIO.read(alignment_path, 'fasta')
    alignments.assert_alignment_length(alignment, length=3)
    alignments.assert_alignment_length(alignment, min=2)
    alignments.assert_alignment_length(alignment, max=4)
    with pytest.raises(AssertionError):
        alignments.assert_alignment_length(alignment, length=1)
        alignments.assert_alignment_length(alignment, min=4)
        alignments.assert_alignment_length(alignment, max=2)
