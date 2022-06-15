import re

import pytest
from Bio.Seq import Seq

from phytest import Sequence
from phytest.utils import PhytestAssertion, PhytestWarning


def test_assert_valid_alphabet():
    sequence = Sequence(
        Seq("ACGTACGTACGT"),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    protein = Sequence(
        Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
        id="PROTEINID",
        name="TEST Protein",
        description="Test protein sequence",
    )
    sequence.assert_valid_alphabet()
    with pytest.raises(
        PhytestAssertion,
        match="Invalid pattern found in 'DNAID'.\nCharacter 'G' at position 3 found which is not in alphabet 'ABCDE'.",
    ):
        sequence.assert_valid_alphabet(alphabet="ABCDE")

    protein.assert_valid_alphabet(alphabet="ACDEFGHIKLMNPQRSTVWYXBZJ")
    with pytest.raises(
        PhytestAssertion,
        match="Invalid pattern found in 'PROTEINID'.\nCharacter 'M' at position 1 found which is not in alphabet 'ATCGN-'.",
    ):
        protein.assert_valid_alphabet()


def test_assert_length():
    sequence = Sequence(
        Seq("A" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_length(length=100, min=99, max=101)

    with pytest.raises(
        PhytestAssertion, match=re.escape("Sequence length of 'DNAID' (100) is not equal to the required length of 1.")
    ):
        sequence.assert_length(length=1)

    with pytest.raises(
        PhytestAssertion, match=re.escape("Sequence length of 'DNAID' (100) is less than the minimum 101.")
    ):
        sequence.assert_length(min=101)

    with pytest.raises(
        PhytestAssertion, match=re.escape("Sequence length of 'DNAID' (100) is greater than the maximum 99.")
    ):
        sequence.assert_length(max=99)

    with pytest.warns(
        PhytestWarning, match=re.escape("Sequence length of 'DNAID' (100) is greater than the maximum 99.")
    ):
        sequence.warn_length(max=99)


def test_assert_count():
    sequence = Sequence(
        Seq("ATG" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_count(pattern='A', count=100, min=99, max=101)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "Sequence 'DNAID' matches pattern 'A' 100 time(s).\nThis is not equal to the required number of 1."
        ),
    ):
        sequence.assert_count(pattern='A', count=1)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape("Sequence 'DNAID' matches pattern 'A' 100 time(s).\nThis is less than the minimum 101."),
    ):
        sequence.assert_count(pattern='A', min=101)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape("Sequence 'DNAID' matches pattern 'A' 100 time(s).\nThis is greater than the maximum 99."),
    ):
        sequence.assert_count(pattern='A', max=99)


def test_assert_count_Ns():
    sequence = Sequence(
        Seq("ATGN" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_count_Ns(count=100, min=99, max=101)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "Sequence 'DNAID' matches pattern 'N' 100 time(s).\nThis is not equal to the required number of 1."
        ),
    ):
        sequence.assert_count_Ns(count=1)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape("Sequence 'DNAID' matches pattern 'N' 100 time(s).\nThis is less than the minimum 101."),
    ):
        sequence.assert_count_Ns(min=101)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape("Sequence 'DNAID' matches pattern 'N' 100 time(s).\nThis is greater than the maximum 99."),
    ):
        sequence.assert_count_Ns(max=99)


def test_assert_count_gaps():
    sequence = Sequence(
        Seq("ATG-" * 100),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_count_gaps(count=100, min=99, max=101)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "Sequence 'DNAID' matches pattern '-' 100 time(s).\nThis is not equal to the required number of 1."
        ),
    ):
        sequence.assert_count_gaps(count=1)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape("Sequence 'DNAID' matches pattern '-' 100 time(s).\nThis is less than the minimum 101."),
    ):
        sequence.assert_count_gaps(min=101)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape("Sequence 'DNAID' matches pattern '-' 100 time(s).\nThis is greater than the maximum 99."),
    ):
        sequence.assert_count_gaps(max=99)


def test_assert_sequence_longest_stretch():
    sequence = Sequence(
        Seq("A" * 10 + "-" * 3 + "N" * 10),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_longest_stretch(pattern='A', count=10, min=9, max=11)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "The longest stretch of pattern 'A' in sequence 'DNAID' is 10.\nThis is not equal to the required number of 1"
        ),
    ):
        sequence.assert_longest_stretch(pattern='A', count=1)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "The longest stretch of pattern 'A' in sequence 'DNAID' is 10.\nThis is less than the minimum 11."
        ),
    ):
        sequence.assert_longest_stretch(pattern='A', min=11)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "The longest stretch of pattern 'A' in sequence 'DNAID' is 10.\nThis is greater than the maximum 9."
        ),
    ):
        sequence.assert_longest_stretch(pattern='A', max=9)


def test_assert_sequence_longest_Ns():
    sequence = Sequence(
        Seq("A" * 10 + "-" * 3 + "N" * 10),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_longest_stretch_Ns(count=10, min=9, max=11)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "The longest stretch of pattern 'N' in sequence 'DNAID' is 10.\nThis is not equal to the required number of 1."
        ),
    ):
        sequence.assert_longest_stretch_Ns(count=1)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "The longest stretch of pattern 'N' in sequence 'DNAID' is 10.\nThis is less than the minimum 11."
        ),
    ):
        sequence.assert_longest_stretch_Ns(min=11)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "The longest stretch of pattern 'N' in sequence 'DNAID' is 10.\nThis is greater than the maximum 9."
        ),
    ):
        sequence.assert_longest_stretch_Ns(max=9)


def test_assert_sequence_longest_gaps():
    sequence = Sequence(
        Seq("A" * 10 + "-" * 3 + "N" * 10),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_longest_stretch_gaps(count=3, min=2, max=4)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "The longest stretch of pattern '-' in sequence 'DNAID' is 3.\nThis is not equal to the required number of 1."
        ),
    ):
        sequence.assert_longest_stretch_gaps(count=1)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "The longest stretch of pattern '-' in sequence 'DNAID' is 3.\nThis is less than the minimum 4."
        ),
    ):
        sequence.assert_longest_stretch_gaps(min=4)
    with pytest.raises(
        PhytestAssertion,
        match=re.escape(
            "The longest stretch of pattern '-' in sequence 'DNAID' is 3.\nThis is greater than the maximum 2."
        ),
    ):
        sequence.assert_longest_stretch_gaps(max=2)


def test_assert_sequence_startswith():
    sequence = Sequence(
        Seq("ATG" + "-" * 3 + "UGA"),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_startswith(pattern='ATG')
    with pytest.raises(PhytestAssertion, match=re.escape("Sequence 'DNAID' does not start with 'UGA'.")):
        sequence.assert_startswith(pattern='UGA')


def test_assert_sequence_endswith():
    sequence = Sequence(
        Seq("ATG" + "-" * 3 + "UGA"),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_endswith(pattern='UGA')
    with pytest.raises(PhytestAssertion, match=re.escape("Sequence 'DNAID' does not end with 'ATG'.")):
        sequence.assert_endswith(pattern='ATG')


def test_assert_sequence_contains():
    sequence = Sequence(
        Seq("ATG" + "TGACGT" + "UGA"),
        id="DNAID",
        name="TEST",
        description="Test dna sequence",
    )
    sequence.assert_contains(pattern='TGACGT')
    with pytest.raises(
        PhytestAssertion,
        match=re.escape("Sequence 'DNAID' matches pattern 'CAGCTG' 0 time(s).\nThis is less than the minimum 1."),
    ):
        sequence.assert_contains(pattern='CAGCTG')
