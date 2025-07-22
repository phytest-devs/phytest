from typing import Iterator

from phytest import Alignment, Data, Sequence, Tree


def test_length(sequence: Sequence):
    sequence.assert_length(length=100)


def test_alignment_length(alignment: Alignment):
    alignment.assert_length(length=4)


def test_tree_number_of_tips(tree: Tree):
    tree.assert_number_of_tips(4)


def test_number_of_trees(trees: Iterator[Tree]):
    assert len(list(trees)) == 1, "There should be exactly one tree in the fixture"


def test_data_number_of_rows(data: Data):
    data.assert_match('name', 'Sequence_[A-D]')
