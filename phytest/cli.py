from pathlib import Path
from typing import Optional

import pytest
import typer

app = typer.Typer()


@app.command()
def cli(
    testfile: Path = typer.Argument(..., help="Path to test file"),
    alignment: Optional[Path] = typer.Option(
        None, "--alignment", "-a", dir_okay=False, exists=True, help="Path to sequence alignment in fasta format"
    ),
    tree: Optional[Path] = typer.Option(None, "--tree", "-t", dir_okay=False, exists=True, help="Path to tree file"),
    report: Optional[bool] = typer.Option(False, "--report", "-r", help="Generate an html report"),
    verbose: Optional[bool] = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    args = [testfile]
    if not verbose:
        args.extend(["-ra", "--tb=no", "--no-header"])
    if alignment is not None:
        args.extend(["--alignment", alignment])
    if tree is not None:
        args.extend(["--tree", tree])
    if report:
        args.extend(["--html=report.html", "--self-contained-html"])

    pytest.main(args, plugins=['phytest'])
