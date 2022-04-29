from io import StringIO

import pytest

from phytest import Tree


def test_assert_tree_number_of_tips():
    treedata = "(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460);"
    handle = StringIO(treedata)
    tree = Tree.read(handle, "newick")
    tree.assert_number_of_tips(tips=7, min=6, max=8)
    with pytest.raises(AssertionError):
        tree.assert_number_of_tips(tips=1)
    with pytest.raises(AssertionError):
        tree.assert_number_of_tips(min=8)
    with pytest.raises(AssertionError):
        tree.assert_number_of_tips(max=6)


def test_assert_tree_is_bifurcating():
    treedata = "(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460);"
    handle = StringIO(treedata)
    tree = Tree.read(handle, "newick")
    tree.assert_is_bifurcating()


def test_assert_tree_total_branch_length():
    treedata = "(Bovine:1,(Hylobates:1,(Pongo:1,(G._Gorilla:1, (P._paniscus:1,H._sapiens:1):1):1):1):1, Rodent:1);"
    handle = StringIO(treedata)
    tree = Tree.read(handle, "newick")
    tree.assert_total_branch_length(length=11, min=10, max=12)
    with pytest.raises(AssertionError):
        tree.assert_total_branch_length(length=1)
    with pytest.raises(AssertionError):
        tree.assert_total_branch_length(min=12)
    with pytest.raises(AssertionError):
        tree.assert_total_branch_length(max=10)
