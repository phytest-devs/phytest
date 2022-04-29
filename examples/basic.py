from phytest import Alignment, Sequence


def test_length(sequence: Sequence):
    sequence.assert_length(sequence, length=100)


def test_alignment_length(alignment: Alignment):
    alignment.assert_length(alignment, length=3)
