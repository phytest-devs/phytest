from phytest import Alignment, Sequence, Tree, bio


def test_sequence(sequence: Sequence):
    bio.sequence.assert_sequence_valid_alphabet(sequence)


def test_alignment(alignment: Alignment):
    bio.alignment.assert_alignment_length(alignment, 52)
    bio.alignment.assert_alignment_width(alignment, 462)


def test_tree(tree: Tree):
    bio.tree.assert_tree_is_bifurcating(tree)
    bio.tree.assert_tree_number_of_tips(tree, 52)
