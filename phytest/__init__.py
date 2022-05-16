from pathlib import Path

import pandas as pd
import pytest
from pandas import DataFrame

from .bio import Alignment, Sequence, Tree
from .main import main as main


def pytest_addoption(parser):
    parser.addoption("--alignment", "-A", action="store", default=None, help="alignment file")
    parser.addoption("--alignment-format", action="store", default='fasta', help="alignment file format")
    parser.addoption("--tree", "-T", action="store", default=None, help="tree file")
    parser.addoption("--tree-format", action="store", default='newick', help="tree file format")
    parser.addoption("--data", "-D", action="store", default=None, help="data file")
    parser.addoption(
        "--apply-fixes", action="store_true", default=None, help="automatically apply fixes where possible"
    )


def pytest_generate_tests(metafunc):
    alignment_path = metafunc.config.getoption("alignment")
    if 'alignment' in metafunc.fixturenames:
        if alignment_path is None:
            raise ValueError(f"{metafunc.function.__name__} requires an alignment file")
        fpth = Path(alignment_path)
        if not fpth.exists():
            raise FileNotFoundError(f"Unable to locate requested alignment file ({fpth})! 😱")
    tree_path = metafunc.config.getoption("tree")
    if 'tree' in metafunc.fixturenames:
        if tree_path is None:
            raise ValueError(f"{metafunc.function.__name__} requires a tree file")
        fpth = Path(tree_path)
        if not fpth.exists():
            raise FileNotFoundError(f"Unable to locate requested tree file ({fpth})! 😱")
        tree_format = metafunc.config.getoption("--tree-format")
        trees = Tree.parse(tree_path, tree_format)
        metafunc.parametrize("tree", trees, ids=lambda t: t.name)
    data_path = metafunc.config.getoption("data")
    if 'data' in metafunc.fixturenames:
        if data_path is None:
            raise ValueError(f"{metafunc.function.__name__} requires a data file")
        fpth = Path(data_path)
        if not fpth.exists():
            raise FileNotFoundError(f"Unable to locate requested data file ({fpth})! 😱")
    if "sequence" in metafunc.fixturenames:
        if alignment_path is None:
            raise ValueError(f"{metafunc.function.__name__} requires an alignment file")
        fpth = Path(alignment_path)
        if not fpth.exists():
            raise FileNotFoundError(f"Unable to locate requested alignment file ({fpth})! 😱")
        alignment_format = metafunc.config.getoption("--alignment-format")
        sequences = Sequence.parse(alignment_path, alignment_format)
        metafunc.parametrize("sequence", sequences, ids=lambda s: s.id)


@pytest.fixture(scope="session", name="alignment")
def _alignment_fixture(request):
    alignment_path = request.config.getoption("alignment")
    alignment_format = request.config.getoption("--alignment-format")
    alignment = Alignment.read(alignment_path, alignment_format)
    return alignment


@pytest.fixture(scope="session", name="data")
def _data_fixture(request):
    data_path = request.config.getoption("data")
    data = pd.read_csv(data_path)
    return data


@pytest.fixture()
def should_fix(request):
    if request.session.testsfailed and request.config.getoption("--apply-fixes"):
        return True
    return False


def pytest_html_report_title(report):
    report.title = "Quality control checks"
