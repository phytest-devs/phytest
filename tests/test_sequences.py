import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phytest import sequences


def test_assert_valid_alphabet():
    dna = SeqRecord(
        Seq("ACGTACGTACGT"),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    protein = SeqRecord(
        Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
        id="PROTEINID",
        name="TEST Protein",
        description="Test protein sequence",
    )
    sequences.assert_sequence_valid_alphabet(dna)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_valid_alphabet(dna, alphabet="ABCDE")

    sequences.assert_sequence_valid_alphabet(protein, alphabet="ACDEFGHIKLMNPQRSTVWYXBZJ")
    with pytest.raises(AssertionError):
        sequences.assert_sequence_valid_alphabet(protein)


def test_assert_length():
    dna = SeqRecord(
        Seq("A" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_length(dna, length=100)
    sequences.assert_sequence_length(dna, min=99)
    sequences.assert_sequence_length(dna, max=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_length(dna, length=1)
        sequences.assert_sequence_length(dna, min=101)
        sequences.assert_sequence_length(dna, max=99)


def test_assert_count():
    dna = SeqRecord(
        Seq("ATG" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_count(dna, base='A', count=100)
    sequences.assert_sequence_count(dna, base='A', min=99)
    sequences.assert_sequence_count(dna, base='A', max=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count(dna, base='A', count=1)
        sequences.assert_sequence_count(dna, base='A', min=101)
        sequences.assert_sequence_count(dna, base='A', max=99)


def test_assert_count_Ns():
    dna = SeqRecord(
        Seq("ATGN" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_count_Ns(dna, count=100)
    sequences.assert_sequence_count_Ns(dna, min=99)
    sequences.assert_sequence_count_Ns(dna, max=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count_Ns(dna, count=1)
        sequences.assert_sequence_count_Ns(dna, min=101)
        sequences.assert_sequence_count_Ns(dna, max=99)


def test_assert_count_gaps():
    dna = SeqRecord(
        Seq("ATG-" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_count_gaps(dna, count=100)
    sequences.assert_sequence_count_gaps(dna, min=99)
    sequences.assert_sequence_count_gaps(dna, max=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count_gaps(dna, count=1)
        sequences.assert_sequence_count_gaps(dna, min=101)
        sequences.assert_sequence_count_gaps(dna, max=99)
