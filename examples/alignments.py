from phytest.alignments import assert_alignment_length
from phytest.sequences import assert_sequence_length


def test_length(sequence):
    assert_sequence_length(sequence, length=100)


def test_alignment_length(alignment):
    assert_alignment_length(alignment, length=3)
