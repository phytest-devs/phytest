import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phytest import alignments

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
    alignments.assert_valid_alphabet(dna)
    with pytest.raises(AssertionError):
        alignments.assert_valid_alphabet(dna, alphabet="ABCDE")    

    alignments.assert_valid_alphabet(protein, alphabet="ACDEFGHIKLMNPQRSTVWYXBZJ")
    with pytest.raises(AssertionError):
        alignments.assert_valid_alphabet(protein)    
