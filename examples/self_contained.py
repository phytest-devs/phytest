import sys

import phytest


def test_length(sequence: phytest.Sequence):
    sequence.assert_length(length=462)


if __name__ == "__main__":
    sys.exit(phytest.main(alignment='examples/data/ice_viruses.fasta', tree='examples/data/ice_viruses.fasta.treefile'))
