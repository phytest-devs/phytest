from phytest import Alignment, Data, Sequence, Tree


def test_alignment_length(alignment: Alignment):
    alignment.assert_length(length=4)
