from pathlib import Path
from typing import Optional

import typer

from .main import main

app = typer.Typer()


@app.command()
def cli(
    testfile: Path = typer.Argument(..., help="Path to test file"),
    alignment: Optional[Path] = typer.Option(
        None, "--alignment", "-a", dir_okay=False, exists=True, help="Path to alignment file"
    ),
    alignment_format: Optional[str] = typer.Option(
        'fasta', "--alignment-format", dir_okay=False, exists=True, help="Alignment format"
    ),
    tree: Optional[Path] = typer.Option(None, "--tree", "-t", dir_okay=False, exists=True, help="Path to tree file"),
    tree_format: Optional[str] = typer.Option(
        'newick', "--tree-format", dir_okay=False, exists=True, help="Tree file format"
    ),
    data: Optional[Path] = typer.Option(None, "--data", "-d", dir_okay=False, exists=True, help="Metadata file (CSV)"),
    report: Optional[Path] = typer.Option(
        None, "--report", "-r", dir_okay=False, exists=False, help="Path to HTML report to generate."
    ),
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
        data=data,
        verbose=verbose,
        report=report,
        expression=expression,
    )
    raise typer.Exit(code=exit_code)
