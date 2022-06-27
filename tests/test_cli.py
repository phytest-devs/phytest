from pathlib import Path

import pytest
from typer.testing import CliRunner

from phytest.cli import app

runner = CliRunner()


def test_cli_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "TESTFILE  Path to test file  [required]" in result.stdout


def test_cli_no_input_file(request: pytest.FixtureRequest):
    result = runner.invoke(app, [str(request.path.parent / "input/testfile1.py")])
    assert result.exit_code == 0
    assert "testfile1.py" in result.stdout
    assert "1 passed" in result.stdout


def test_cli_basic(request: pytest.FixtureRequest):
    result = runner.invoke(
        app,
        [
            str(request.path.parent / "input/basic.py"),
            "-s",
            "examples/data/example.fasta",
            "-t",
            "examples/data/example.tree",
            "-d",
            "examples/data/example.csv",
        ],
    )
    assert "7 passed" in result.stdout


def test_cli_report(request: pytest.FixtureRequest):
    result = runner.invoke(
        app,
        [
            str(request.path.parent / "input/basic.py"),
            "-s",
            "examples/data/example.fasta",
            "-t",
            "examples/data/example.tree",
            "-d",
            "examples/data/example.csv",
            "-r",
            "pytest-report.html",
        ],
    )
    assert Path("pytest-report.html").exists()
