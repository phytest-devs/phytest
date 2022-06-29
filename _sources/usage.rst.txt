==============
Running Tests
==============

Phytest is based on the popular Python testing framework Pytest.
It provides many convenient helper functions for testing phylogenetic analyses including methods for testing sequences, alignments, trees and metadata files.
Phytest has been developed as a command-line interface, Python module, and Pytest plugin, providing multiple methods of invocation.

Phytest CLI
===========

Phytest provides a command line interface (CLI) for running testing on specific data files.

.. code-block:: python

    from phytest import Sequence

    def test_gc_content(sequence: Sequence):
        sequence.assert_percent_GC(
            min=30,
            max=40
        )

The Phytest CLI requires a path to the file containing user defined tests and has optional flags for specifying sequence/alignment, tree and metadata files:

.. code-block:: bash

    phytest test.py --sequence sequences.fasta --tree tree.newick --data metadata.csv

Alternative file formats can be specified with :code:`--sequence-format`, :code:`--tree-format`, :code:`--data-format` flags.
Supported formats include those supported by Biopython (sequences and trees) and TSV and Excel (metadata).

Phytest Module
================

The Phytest module can be imported into script so that tests are self-contained i.e. data files are specified in the tests.

.. code-block:: python

    import phytest

    def test_gc_content(sequence: phytest.Sequence):
        sequence.assert_percent_GC(
            min=30,
            max=40
        )

    if __name__ == "__main__":
        sys.exit(phytest.main(sequence='examples/data/ice_viruses.fasta'))

This test style uses :code:`if __name__ == "__main__"` python convention to only run the tests when invoked from the command line using the python command.

.. code-block:: bash

    python test.py

The :code:`phytest.main` function will run the tests and return a exit status (0 ir 1) that is passed to :code:`sys.exit` to ensure the test exit correctly.


Pytest Plugin
================

Phytest can also be used as a Pytest plugin. Simply install Phytest and then run Pytest on the test file with the appropriate flags.

.. code-block:: bash

    pytest test.py --sequence sequences.fasta

.. NOTE::
   Short hand flags must be capitalised when running Phytest through Pytest e.g. :code:`pytest test.py -S sequences.fasta`
