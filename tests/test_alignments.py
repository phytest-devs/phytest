import pytest

from phytest import Alignment


def test_assert_alignment_width():
    alignment_path = 'examples/data/invalid.fasta'
    alignment = Alignment.read(alignment_path, 'fasta')
    alignment.assert_width(width=100, min=99, max=101)
    with pytest.raises(AssertionError):
        alignment.assert_width(width=99)
    with pytest.raises(AssertionError):
        alignment.assert_width(min=101)
    with pytest.raises(AssertionError):
        alignment.assert_width(max=99)


def test_assert_alignment_length():
    alignment_path = 'examples/data/invalid.fasta'
    alignment = Alignment.read(alignment_path, 'fasta')
    alignment.assert_length(length=3, min=2, max=4)
    with pytest.raises(AssertionError):
        alignment.assert_length(length=1)
    with pytest.raises(AssertionError):
        alignment.assert_length(min=4)
    with pytest.raises(AssertionError):
        alignment.assert_length(max=2)
