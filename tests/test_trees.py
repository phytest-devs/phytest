from io import StringIO

import pytest
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phytest import tree


def test_assert_tree_number_of_tips():
    treedata = "(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460);"
    handle = StringIO(treedata)
    tree = Phylo.read(handle, "newick")
    tree.assert_tree_number_of_tips(tree, tips=7, min=6, max=8)
    with pytest.raises(AssertionError):
        tree.assert_tree_number_of_tips(tree, tips=1)
    with pytest.raises(AssertionError):
        tree.assert_tree_number_of_tips(tree, min=8)
    with pytest.raises(AssertionError):
        tree.assert_tree_number_of_tips(tree, max=6)


def test_assert_tree_is_bifurcating():
    treedata = "(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460);"
    handle = StringIO(treedata)
    tree = Phylo.read(handle, "newick")
    tree.assert_tree_is_bifurcating(tree)


def test_assert_tree_total_branch_length():
    treedata = "(Bovine:1,(Hylobates:1,(Pongo:1,(G._Gorilla:1, (P._paniscus:1,H._sapiens:1):1):1):1):1, Rodent:1);"
    handle = StringIO(treedata)
    tree = Phylo.read(handle, "newick")
    tree.assert_tree_total_branch_length(tree, length=11, min=10, max=12)
    with pytest.raises(AssertionError):
        tree.assert_tree_total_branch_length(tree, length=1)
    with pytest.raises(AssertionError):
        tree.assert_tree_total_branch_length(tree, min=12)
    with pytest.raises(AssertionError):
        tree.assert_tree_total_branch_length(tree, max=10)
