from pathlib import Path
from typing import Optional

import typer

from .main import main

app = typer.Typer()


@app.command()
def cli(
    testfile: Path = typer.Argument(..., help="Path to test file"),
    sequence: Optional[Path] = typer.Option(
        None, "--sequence", "-s", dir_okay=False, exists=True, help="Path to sequence file"
    ),
    sequence_format: Optional[str] = typer.Option(
        'fasta', "--sequence-format", dir_okay=False, exists=True, help="Sequence file format"
    ),
    tree: Optional[Path] = typer.Option(None, "--tree", "-t", dir_okay=False, exists=True, help="Path to tree file"),
    tree_format: Optional[str] = typer.Option(
        'newick', "--tree-format", dir_okay=False, exists=True, help="Tree file format"
    ),
    data: Optional[Path] = typer.Option(None, "--data", "-d", dir_okay=False, exists=True, help="Metadata file (CSV)"),
    data_format: Optional[str] = typer.Option(
        'csv', "--data-format", dir_okay=False, exists=True, help="Data file format"
    ),
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
        sequence=sequence,
        sequence_format=sequence_format,
        tree=tree,
        tree_format=tree_format,
        data=data,
        data_format=data_format,
        verbose=verbose,
        report=report,
        expression=expression,
    )
    raise typer.Exit(code=exit_code)
