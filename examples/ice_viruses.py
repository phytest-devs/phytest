from phytest import alignments, sequences, trees


def test_sequence(sequence):
    sequences.assert_sequence_valid_alphabet(sequence)


def test_alignment(alignment):
    alignments.assert_alignment_length(alignment, 52)
    alignments.assert_alignment_width(alignment, 462)


def test_tree(tree):
    trees.assert_tree_is_bifurcating(tree)
    trees.assert_tree_number_of_tips(tree, 52)
