from itertools import zip_longest
from pathlib import Path

import pytest


def pytest_addoption(parser):
    parser.addoption("--alignments", "-A", action="store", default=None, help="alignments fasta file")
    parser.addoption("--tree", "-T", action="store", default=None, help="tree file")
    parser.addoption(
        "--apply-fixes", action="store_true", default=None, help="automatically apply fixes where possible"
    )


def load_sequences(fpth: Path):
    with fpth.open("r") as fp:
        for seq_id, seq in zip_longest(fp, fp):
            yield seq_id.strip(), seq.strip()


def pytest_generate_tests(metafunc):
    if "sequence" in metafunc.fixturenames:
        alignments_path = metafunc.config.getoption("alignments")
        if alignments_path is None:
            pytest.skip(f"{metafunc.function.__name__} requires alignments file")
        fpth = Path(alignments_path)
        if not fpth.exists():
            raise FileNotFoundError("Unable to locate requested input file! 😱")
        sequences = load_sequences(fpth)
        metafunc.parametrize("sequence", sequences, ids=lambda s: s[0].lstrip('>'))


@pytest.fixture()
def should_fix(request):
    if request.session.testsfailed and request.config.getoption("--apply-fixes"):
        return True
    return False


def pytest_html_report_title(report):
    report.title = "Quality control checks"
