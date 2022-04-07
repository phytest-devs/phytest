![POC meme](https://memegenerator.net/img/instances/67223136.jpg)

# PhyPyTest POC

## General info and questions

This is an example of how we could use pytest for the quality control checks.

The idea would be to take `./phypytest/alignments/*` as an example and then wrap some of that functionality up with utility functions to make it reusable. We could then wrap the whole thing in a CLI and provide a library of tests for sequences, alignments and trees etc.

## Install

```
poetry install
```

## Current example (no custom CLI yet)

```
poetry run pytest --tb=line phypytest/alignments/test_alignments.py --input examples/data/invalid.fasta -s --apply-fixes
```

Generate a report by adding:

```
--html=report.html --self-contained-html
```

this report can be customised in future (see the [pytest-html user guide](https://pytest-html.readthedocs.io/en/latest/user_guide.html)).


## Still to do

_...even if we are happy with the interface._

- [ ] Better stdout reporting if we can.
- [ ] Test groups for keeping related tests together in report.
- [ ] Wrap up common functionality into utility functions etc.
- [ ] A suite of builtin tests for sequencing and tree data
- [ ] Improve html report output and add metadata
- [ ] a million other things


## Other considerations

This is more complicated than [the phytest POC](https://gitlab.unimelb.edu.au/mdap/phytest), but way more feature full. We can also run tests in parallel etc. without worrying about the state hacks in phytest. In reality maybe we will want somewhere between the two. Or perhaps to go down the phytest path but with a less hacky interface that uses class inheritance. 🤷
