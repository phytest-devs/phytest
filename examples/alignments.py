from phytest.alignments import assert_alignment_length, assert_length


def test_length(sequence):
    assert_length(sequence, length=100)


def test_alignment_length(alignment):
    assert_alignment_length(alignment, length=3)
