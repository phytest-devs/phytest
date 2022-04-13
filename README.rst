
==============
phytest
==============

.. start-badges

|pipline badge| |coverage badge| |docs badge| |black badge|

.. |pipline badge| image:: https://github.com/smutch/phytest/actions/workflows/coverage.yml/badge.svg
    :target: https://github.com/smutch/phytest/actions

.. |docs badge| image:: https://github.com/smutch/phytest/actions/workflows/docs.yml/badge.svg
    :target: https://smutch.github.io/phytest/
    
.. |black badge| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    
.. |coverage badge| image:: https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/smutch/506563cd9b49c8126284e34864c862d0/raw/coverage-badge.json
    :target: https://smutch.github.io/phytest/coverage/
    
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

.. code-block:: bash

    phytest phytest/alignments/test_alignments.py -A examples/data/invalid.fasta

Generate a report by adding:

.. code-block:: bash

    --html=report.html --self-contained-html

This report can be customised in future (see the `pytest-html user guide <https://pytest-html.readthedocs.io/en/latest/user_guide.html>`_.

More information to come.

.. end-quickstart


Old info
========


This is an example of how we could use pytest for the quality control checks.

The idea would be to take `./phytest/alignments/*` as an example and then wrap some of that functionality up with utility functions to make it reusable. We could then wrap the whole thing in a CLI and provide a library of tests for sequences, alignments and trees etc.


Still to do
====================

_...even if we are happy with the interface._

- [ ] Better stdout reporting if we can.
- [ ] Test groups for keeping related tests together in report.
- [ ] Wrap up common functionality into utility functions etc.
- [ ] A suite of builtin tests for sequencing and tree data
- [ ] Improve html report output and add metadata
- [ ] a million other things


Other considerations
====================

This is more complicated than [the phytest POC](https://gitlab.unimelb.edu.au/mdap/phytest), but way more feature full. We can also run tests in parallel etc. without worrying about the state hacks in phytest. In reality maybe we will want somewhere between the two. Or perhaps to go down the phytest path but with a less hacky interface that uses class inheritance. 🤷