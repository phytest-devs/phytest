from phytest import Alignment, Sequence, asserts


def test_length(sequence: Sequence):
    asserts.sequences.assert_sequence_length(sequence, length=100)


def test_alignment_length(alignment: Alignment):
    asserts.alignments.assert_alignment_length(alignment, length=3)
