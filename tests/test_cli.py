from typer.testing import CliRunner

from phytest.cli import app

runner = CliRunner()


def test_cli_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "TESTFILE  Path to test file  [required]" in result.stdout
