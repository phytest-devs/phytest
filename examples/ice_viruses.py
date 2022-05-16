from phytest import Alignment, Sequence, Tree


def test_sequence(sequence: Sequence):
    sequence.assert_valid_alphabet()


def test_alignment(alignment: Alignment):
    alignment.assert_width(462)


def test_root_to_tip(tree: Tree):
    # Evolutionary rate of influenza is on the order of 10^-3
    tree.assert_root_to_tip(min_r_squared=0.5, min_rate=0.001, max_rate=0.009, min_root_date=1852, max_root_date=1853)
