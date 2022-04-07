from itertools import zip_longest
from pathlib import Path


def pytest_addoption(parser):
    parser.addoption("--input", action="store", default=None, help="input fasta file")


def load_sequences(fpth: Path):
    with fpth.open("r") as fp:
        for seq_id, seq in zip_longest(fp, fp):
            yield seq_id.strip(), seq.strip()


def pytest_generate_tests(metafunc):
    if "sequence" in metafunc.fixturenames:
        fpth = Path(metafunc.config.getoption("input"))
        print(fpth)
        if not fpth.exists():
            raise FileNotFoundError("Unable to locate requested input file! 😱")
        sequences = load_sequences(fpth)
        metafunc.parametrize("sequence", sequences, ids=lambda s: s[0].lstrip('>'))


def pytest_html_report_title(report):
    report.title = "Quality control checks"
