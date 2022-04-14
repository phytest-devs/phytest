from pathlib import Path
from typing import Optional

import pytest
import typer

app = typer.Typer()


@app.command()
def cli(
    testfile: Path = typer.Argument(..., help="Path to test file"),
    alignment: Optional[Path] = typer.Option(
        None, "--alignment", "-a", dir_okay=False, exists=True, help="Path to sequence alignment"
    ),
    alignment_format: Optional[str] = typer.Option(
        'fasta', "--alignment-format", dir_okay=False, exists=True, help="Sequence alignment format"
    ),
    tree: Optional[Path] = typer.Option(None, "--tree", "-t", dir_okay=False, exists=True, help="Tree file format"),
    tree_format: Optional[str] = typer.Option(
        'newick', "--tree-format", dir_okay=False, exists=True, help="Sequence alignment format"
    ),
    report: Optional[bool] = typer.Option(False, "--report", "-r", help="Generate an html report"),
    verbose: Optional[bool] = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    args = [testfile]
    if not verbose:
        args.extend(["-ra", "--tb=no", "--no-header"])
    if alignment is not None:
        args.extend(["--alignment", alignment])
        args.extend(["--alignment-format", alignment_format])
    if tree is not None:
        args.extend(["--tree", tree])
        args.extend(["--tree-format", tree_format])
    if report:
        args.extend(["--html=report.html", "--self-contained-html"])

    pytest.main(args, plugins=['phytest'])
