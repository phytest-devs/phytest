import phytest
from phytest import sequences


def test_length(sequence):
    sequences.assert_sequence_length(sequence, length=462)


if __name__ == "__main__":
    phytest.main(alignment='examples/data/ice_viruses.fasta', tree='examples/data/ice_viruses.fasta.treefile')
