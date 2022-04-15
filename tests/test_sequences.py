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
    sequences.assert_sequence_length(dna, length=100, min=99, max=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_length(dna, length=1)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_length(dna, min=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_length(dna, max=99)


def test_assert_count():
    dna = SeqRecord(
        Seq("ATG" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_count(dna, pattern='A', count=100, min=99, max=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count(dna, pattern='A', count=1)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count(dna, pattern='A', min=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count(dna, pattern='A', max=99)


def test_assert_count_Ns():
    dna = SeqRecord(
        Seq("ATGN" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_count_Ns(dna, count=100, min=99, max=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count_Ns(dna, count=1)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count_Ns(dna, min=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count_Ns(dna, max=99)


def test_assert_count_gaps():
    dna = SeqRecord(
        Seq("ATG-" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_count_gaps(dna, count=100, min=99, max=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count_gaps(dna, count=1)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count_gaps(dna, min=101)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_count_gaps(dna, max=99)


def test_assert_sequence_longest_stretch():
    dna = SeqRecord(
        Seq("A" * 10 + "-" * 3 + "N" * 10),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_longest_stretch(dna, pattern='A', count=10, min=9, max=11)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_longest_stretch(dna, pattern='A', count=1)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_longest_stretch(dna, pattern='A', min=11)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_longest_stretch(dna, pattern='A', max=9)


def test_assert_sequence_longest_Ns():
    dna = SeqRecord(
        Seq("A" * 10 + "-" * 3 + "N" * 10),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_longest_stretch_Ns(dna, count=10, min=9, max=11)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_longest_stretch_Ns(dna, count=1)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_longest_stretch_Ns(dna, min=11)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_longest_stretch_Ns(dna, max=9)


def test_assert_sequence_longest_gaps():
    dna = SeqRecord(
        Seq("A" * 10 + "-" * 3 + "N" * 10),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_longest_stretch_gaps(dna, count=3, min=2, max=4)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_longest_stretch_gaps(dna, count=1)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_longest_stretch_gaps(dna, min=4)
    with pytest.raises(AssertionError):
        sequences.assert_sequence_longest_stretch_gaps(dna, max=2)


def test_assert_sequence_startswith():
    dna = SeqRecord(
        Seq("ATG" + "-" * 3 + "UGA"),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_startswith(dna, pattern='ATG')
    with pytest.raises(AssertionError):
        sequences.assert_sequence_startswith(dna, pattern='UGA')


def test_assert_sequence_endswith():
    dna = SeqRecord(
        Seq("ATG" + "-" * 3 + "UGA"),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_endswith(dna, pattern='UGA')
    with pytest.raises(AssertionError):
        sequences.assert_sequence_endswith(dna, pattern='ATG')


def test_assert_sequence_contains_motif():
    dna = SeqRecord(
        Seq("ATG" + "TGACGT" + "UGA"),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequences.assert_sequence_contains_motif(dna, motif='TGACGT')
    with pytest.raises(AssertionError):
        sequences.assert_sequence_contains_motif(dna, motif='CAGCTG')
