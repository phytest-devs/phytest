from typing import Optional
from warnings import warn

from Bio import Phylo as Phylo
from Bio.Phylo.BaseTree import Tree as BioTree


class Tree(BioTree):
    @classmethod
    def read(cls, tree_path, tree_format) -> 'Tree':
        tree = Phylo.read(tree_path, tree_format)
        return Tree(root=tree.root, rooted=tree.rooted, id=tree.id, name=tree.name)

    def assert_number_of_tips(
        self,
        tips: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ):
        number_of_tips = len(self.get_terminals())
        if tips is not None:
            assert number_of_tips == tips
        if min is not None:
            assert number_of_tips >= min
        if max is not None:
            assert number_of_tips <= max
        if warning is not None and number_of_tips != warning:
            warn(f"Number of tips '{number_of_tips}' != {warning}")

    def assert_is_bifurcating(self):
        assert self.is_bifurcating()

    def assert_total_branch_length(
        self,
        length: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ):
        total_branch_length = self.total_branch_length()
        if length is not None:
            assert total_branch_length == length
        if min is not None:
            assert total_branch_length >= min
        if max is not None:
            assert total_branch_length <= max
        if warning is not None and total_branch_length != warning:
            warn(f"Total branch length '{total_branch_length}' != {warning}")
