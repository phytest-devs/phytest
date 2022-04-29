from phytest import Alignment, Sequence, Tree


def test_sequence(sequence: Sequence):
    sequence.assert_valid_alphabet()


def test_alignment(alignment: Alignment):
    alignment.assert_length(52)
    alignment.assert_width(462)


def test_tree(tree: Tree):
    tree.assert_is_bifurcating()
    tree.assert_number_of_tips(52)
