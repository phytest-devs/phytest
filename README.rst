
==============
phytest
==============

.. start-badges

|pipeline badge| |coverage badge| |docs badge| |black badge| |pre-commit badge| |ice-virus-example badge|

.. |pipeline badge| image:: https://github.com/smutch/phytest/workflows/pipeline/badge.svg
    :target: https://github.com/smutch/phytest/actions

.. |docs badge| image:: https://github.com/smutch/phytest/workflows/docs/badge.svg
    :target: https://smutch.github.io/phytest/

.. |black badge| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

.. |coverage badge| image:: https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/smutch/e8160655e03d9015b1e93b97ed611f4f/raw/coverage-badge.json
    :target: https://smutch.github.io/phytest/coverage/

.. |pre-commit badge| image:: https://results.pre-commit.ci/badge/github/phytest-devs/phytest/main.svg
    :target: https://results.pre-commit.ci/latest/github/phytest-devs/phytest/main

.. |ice-virus-example badge| image:: https://github.com/smutch/phytest/workflows/ice-virus-example/badge.svg
    :target: https://github.com/smutch/phytest/actions/workflows/ice_virus_example.yml

.. end-badges

Quality control for phylogenetic pipelines using pytest

.. start-quickstart

Installation
============
Install phytest using pip:

.. code-block:: bash

    pip install git+https://github.com/smutch/phytest.git

.. note ::

    Soon installation will be possible using PyPI.


Usage
============

Phytest will run user defined tests against an alignment and tree file. Here we will create example data files to run our tests on.

Create an alignment fasta file `example.fasta`

.. code-block:: text

    >A
    ATGAGATCCCCGATAGCGAGCTAGCGATCGCAGCGACTCAGCAGCTACAGCGCAGAGGAGAGAGAGGCCCCTATTTACTAGAGCTCCAGATATAGNNTAG
    >B
    ATGAGATCCCCGATAGCGAGCTAGXGATCGCAGCGACTCAGCAGCTACAGCGCAGAGGAGAGAGAGGCCCCTATTTACTAGAGCTCCAGATATAGNNTAG
    >C
    ATGAGA--CCCGATAGCGAGCTAGCGATCGCAGCGACTCAGCAGCTACAGCGCAGAGGAGAGAGAGGCCCCTATTTACTAGAGCTCCAGATATAGNNTAG
    >D
    ATGAGATCCCCGATAGCGAGCTAGCGATCGCAGCGACTCAGCAGCTACAGCGCAGAGGAGAGAGAGGCCCCTATTTACTAGAGCTCCAGATATAGNNTAG
    >E
    ATGAGATCCCCGATAGCGAGCTAGCGATNNNNNNNNNNNNNNNNNTACAGCGCAGAGGAGAGAGAGGCCCCTATTTACTAGAGCTCCAGATATAGNNTAG

Create a tree newick file `example.tree`

.. code-block:: text

    (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);

Writing a test file
########

We want to enforce the follow constraints on our data:
    1. The alignment has 4 sequences
    2. The sequences have a length of 100
    3. The sequences only contains the characters A, T, G, C, N and -
    4. The sequences are allowed to only contain single base deletions
    5. The longest stretch of Ns is 10
    6. The tree has 4 tips
    7. The tree is bifurcating
    8. There are no outlier branches in the tree

We can write these tests in a python files `example.py`

.. code-block:: python

    from phytest import Alignment, Sequence, Tree, asserts
    

    def test_alignment_has_4_sequences(alignment: Alignment):
        asserts.alignments.assert_alignment_length(alignment, 4)


    def test_alignment_has_a_width_of_100(alignment: Alignment):
        asserts.alignments.assert_alignment_width(alignment, 100)


    def test_sequences_only_contains_the_characters(sequence: Sequence):
        asserts.sequences.assert_sequence_valid_alphabet(sequence, alphabet="ATGCN-")


    def test_single_base_deletions(sequence: Sequence):
        asserts.sequences.assert_sequence_longest_stretch_gaps(sequence, max=1)


    def test_longest_stretch_of_Ns_is_10(sequence: Sequence):
        asserts.sequences.assert_sequence_longest_stretch_Ns(sequence, max=10)


    def test_tree_has_4_tips(tree: Tree):
        asserts.trees.assert_tree_number_of_tips(tree, 4)


    def test_tree_is_bifurcating(tree: Tree):
        asserts.trees.assert_tree_is_bifurcating(tree)


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


We can then run these test on our data with `phytest`:

.. code-block:: bash

    phytest examples/example.py -a examples/data/example.fasta -t examples/data/example.tree

Generate a report by adding `--report`.

.. image:: docs/images/report.png
  :alt: HTML Report

This report can be customised in future (see the `pytest-html user guide <https://pytest-html.readthedocs.io/en/latest/user_guide.html>`_).

From the output we can see several tests failed:

.. code-block:: bash

    FAILED examples/example.py::test_sequences_only_contains_the_characters[B] - AssertionError: Invalid pattern found in 'B'!
    FAILED examples/example.py::test_single_base_deletions[C] - AssertionError: Longest stretch of '-' in 'C' > 1!
    FAILED examples/example.py::test_longest_stretch_of_Ns_is_10[D] - AssertionError: Longest stretch of 'N' in 'D' > 10!
    FAILED examples/example.py::test_no_outlier_branches - AssertionError: Outlier tip 'A' (branch length = 1.0)!

    Results (0.07s):
        30 passed
        4 failed
            - examples/example.py:20 test_sequences_only_contains_the_characters[B]
            - examples/example.py:23 test_single_base_deletions[C]
            - examples/example.py:26 test_longest_stretch_of_Ns_is_10[D]
            - examples/example.py:35 test_no_outlier_branches

.. end-quickstart
