from phytest import Alignment, Sequence, bio


def test_length(sequence: Sequence):
    bio.sequence.assert_sequence_length(sequence, length=100)


def test_alignment_length(alignment: Alignment):
    bio.alignment.assert_alignment_length(alignment, length=3)
