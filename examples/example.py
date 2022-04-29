# We want to enforce the follow constraints on our data:
#     1. The alignment has 5 sequences
#     2. The alignment has a width of 100
#     3. The sequences only contains the characters A, T, G, C, N and -
#     4. The sequences are allowed to only contain single base deletions
#     5. The longest stretch of Ns is 10
#     6. The tree has 5 tips
#     7. The tree is bifurcating
#     8. There are no outlier branches in the tree

from phytest import Alignment, Sequence, Tree, bio
from phytest.bio import tree


def test_alignment_has_4_sequences(alignment: Alignment):
    bio.alignment.assert_alignment_length(alignment, 4)


def test_alignment_has_a_width_of_100(alignment: Alignment):
    bio.alignment.assert_alignment_width(alignment, 100)


def test_sequences_only_contains_the_characters(sequence: Sequence):
    bio.sequence.assert_sequence_valid_alphabet(sequence, alphabet="ATGCN-")


def test_single_base_deletions(sequence: Sequence):
    bio.sequence.assert_sequence_longest_stretch_gaps(sequence, max=1)


def test_longest_stretch_of_Ns_is_10(sequence: Sequence):
    bio.sequence.assert_sequence_longest_stretch_Ns(sequence, max=10)


def test_tree_has_4_tips(tree: Tree):
    bio.tree.assert_tree_number_of_tips(tree, 4)


def test_tree_is_bifurcating(tree: Tree):
    bio.tree.assert_tree_is_bifurcating(tree)


def test_no_outlier_branches(tree: Tree):
    # Here we create custom functions to detect outliers
    import statistics

    def get_parent(tree, child_clade):
        node_path = tree.get_path(child_clade)
        if len(node_path) == 1:
            return tree.root
        return node_path[-2]

    branch_lengths = [tree.distance(tip, get_parent(tree, tip)) for tip in tree.get_terminals()]
    for branch_length, tip in zip(branch_lengths, tree.get_terminals()):
        assert branch_length < statistics.mean(branch_lengths) + statistics.stdev(
            branch_lengths
        ), f"Outlier tip '{tip.name}' (branch length = {branch_length})!"
