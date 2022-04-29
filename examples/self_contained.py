import phytest
from phytest.bio.sequence import Sequence


def test_length(sequence: Sequence):
    sequence.assert_length(length=462)


if __name__ == "__main__":
    phytest.main(alignment='examples/data/ice_viruses.fasta', tree='examples/data/ice_viruses.fasta.treefile')
