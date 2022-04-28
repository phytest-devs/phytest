from phytest import Alignment, Sequence, alignments, sequences


def test_length(sequence: Sequence):
    sequences.assert_sequence_length(sequence, length=100)


def test_alignment_length(alignment: Alignment):
    alignments.assert_alignment_length(alignment, length=3)
