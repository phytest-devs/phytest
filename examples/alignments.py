from phytest.alignments import assert_length


def test_length(sequence):
    seq_id, seq = sequence
    assert_length(seq, length=100)
