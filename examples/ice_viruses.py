from phytest import Alignment, Sequence, Tree, asserts


def test_sequence(sequence: Sequence):
    asserts.sequences.assert_sequence_valid_alphabet(sequence)


def test_alignment(alignment: Alignment):
    asserts.alignments.assert_alignment_length(alignment, 52)
    asserts.alignments.assert_alignment_width(alignment, 462)


def test_tree(tree: Tree):
    asserts.trees.assert_tree_is_bifurcating(tree)
    asserts.trees.assert_tree_number_of_tips(tree, 52)
