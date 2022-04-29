from phytest import Alignment, Sequence, Tree


def test_alignment_has_4_sequences(alignment: Alignment):
    alignment.assert_length(4)


def test_alignment_has_a_width_of_100(alignment: Alignment):
    alignment.assert_width(100)


def test_sequences_only_contains_the_characters(sequence: Sequence):
    sequence.assert_valid_alphabet(alphabet="ATGCN-")


def test_single_base_deletions(sequence: Sequence):
    sequence.assert_longest_stretch_gaps(max=1)


def test_longest_stretch_of_Ns_is_10(sequence: Sequence):
    sequence.assert_longest_stretch_Ns(max=10)


def test_tree_has_4_tips(tree: Tree):
    tree.assert_number_of_tips(4)


def test_tree_is_bifurcating(tree: Tree):
    tree.assert_is_bifurcating()


def test_outlier_branches(tree: Tree):
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
