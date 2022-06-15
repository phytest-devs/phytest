import re
import pytest

from phytest import Alignment

from phytest.utils import PhytestAssertion, PhytestWarning

def test_assert_alignment_width():
    alignment_path = 'examples/data/invalid.fasta'
    alignment = Alignment.read(alignment_path, 'fasta')
    alignment.assert_width(width=100, min=99, max=101)
    with pytest.raises(
        PhytestAssertion, 
        match=re.escape("The width of the alignment is 100.\nThis is not equal to the required width of 99.")
    ):
        alignment.assert_width(width=99)
    with pytest.raises(
        PhytestAssertion, 
        match=re.escape("The width of the alignment is 100.\nThis is less than the minimum width of 101.")
    ):
        alignment.assert_width(min=101)
    with pytest.raises(
        PhytestAssertion, 
        match=re.escape("The width of the alignment is 100.\nThis is greater than the maximum width of 99.")
    ):
        alignment.assert_width(max=99)

    with pytest.warns(
        PhytestWarning, 
        match=re.escape("The width of the alignment is 100.\nThis is greater than the maximum width of 99.")
    ):
        alignment.assert_width(max=99, warning=True)


def test_assert_alignment_length():
    alignment_path = 'examples/data/invalid.fasta'
    alignment = Alignment.read(alignment_path, 'fasta')
    alignment.assert_length(length=3, min=2, max=4)
    with pytest.raises(
        PhytestAssertion, 
        match=re.escape("The number of sequences in the alignment is 3.\nThis is less than required number of 1.")
    ):
        alignment.assert_length(length=1)
    with pytest.raises(
        PhytestAssertion, 
        match=re.escape("The number of sequences in the alignment is 3.\nThis is less than the minimum 4.")
    ):
        alignment.assert_length(min=4)
    with pytest.raises(
        PhytestAssertion, 
        match=re.escape("The number of sequences in the alignment is 3.\nThis is greater than the maximum 2.")
    ):
        alignment.assert_length(max=2)
