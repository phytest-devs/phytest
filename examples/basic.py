from phytest import alignments, sequences


def test_length(sequence):
    sequences.assert_sequence_length(sequence, length=100)


def test_alignment_length(alignment):
    alignments.assert_alignment_length(alignment, length=3)
