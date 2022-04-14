from typing import Optional
from warnings import warn

from Bio.Phylo.BaseTree import Tree


def assert_tree_number_of_tips(
    tree: Tree,
    tips: Optional[int] = None,
    *,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
):
    number_of_tips = len(tree.get_terminals())
    if tips is not None:
        assert number_of_tips == tips
    if min is not None:
        assert number_of_tips > min
    if max is not None:
        assert number_of_tips < max
    if warning is not None and number_of_tips != warning:
        warn(f"Number of tips '{number_of_tips}' != {warning}")


def assert_tree_is_bifurcating(tree: Tree):
    assert tree.is_bifurcating()


def assert_tree_total_branch_length(
    tree: Tree,
    length: Optional[int] = None,
    *,
    min: Optional[int] = None,
    max: Optional[int] = None,
    warning: Optional[int] = None,
):
    total_branch_length = tree.total_branch_length()
    if length is not None:
        assert total_branch_length == length
    if min is not None:
        assert total_branch_length > min
    if max is not None:
        assert total_branch_length < max
    if warning is not None and total_branch_length != warning:
        warn(f"Total branch length '{total_branch_length}' != {warning}")
