from itertools import zip_longest
from pathlib import Path

import pytest
from Bio import AlignIO, Phylo, SeqIO

from .helpers import alignments as alignments
from .helpers import sequences as sequences


def pytest_addoption(parser):
    parser.addoption("--alignment", "-A", action="store", default=None, help="alignment fasta file")
    parser.addoption("--tree", "-T", action="store", default=None, help="tree file")
    parser.addoption(
        "--apply-fixes", action="store_true", default=None, help="automatically apply fixes where possible"
    )


def pytest_generate_tests(metafunc):
    alignment_path = metafunc.config.getoption("alignment")
    if 'alignment' in metafunc.fixturenames:
        if alignment_path is None:
            raise ValueError(f"{metafunc.function.__name__} requires alignment file")
    tree_path = metafunc.config.getoption("tree")
    if 'tree' in metafunc.fixturenames:
        if tree_path is None:
            raise ValueError(f"{metafunc.function.__name__} requires tree file")
    if "sequence" in metafunc.fixturenames:
        if alignment_path is None:
            raise ValueError(f"{metafunc.function.__name__} requires alignment file")
        fpth = Path(alignment_path)
        if not fpth.exists():
            raise FileNotFoundError("Unable to locate requested input file! 😱")
        sequences = SeqIO.parse(alignment_path, 'fasta')
        metafunc.parametrize("sequence", sequences, ids=lambda s: s.id)


@pytest.fixture(scope="session", name="alignment")
def _alignment_fixture(request):
    alignment_path = request.config.getoption("alignment")
    try:
        alignment = AlignIO.read(alignment_path, 'fasta')
    except TypeError:
        pass
    return alignment


@pytest.fixture(scope="session", name="tree")
def _tree_fixture(request):
    tree_path = request.config.getoption("tree")
    try:
        tree = Phylo.read(tree_path, "newick")
    except TypeError:
        pass
    return tree


@pytest.fixture()
def should_fix(request):
    if request.session.testsfailed and request.config.getoption("--apply-fixes"):
        return True
    return False


def pytest_html_report_title(report):
    report.title = "Quality control checks"
