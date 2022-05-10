from pathlib import Path
from typing import Optional

import typer

from .main import main

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
    expression: Optional[str] = typer.Option(
        None, "-k", help="Only run tests which match the given substring expression."
    ),
):
    exit_code = main(
        testfile=testfile,
        alignment=alignment,
        alignment_format=alignment_format,
        tree=tree,
        tree_format=tree_format,
        verbose=verbose,
        report=report,
        expression=expression,
    )
    raise typer.Exit(code=exit_code)
