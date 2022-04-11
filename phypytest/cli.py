from pathlib import Path
from typing import Optional

import pytest
import typer


def cli(
    testfile: Path = typer.Argument(..., help="Path to test file"),
    alignments: Optional[Path] = typer.Option(
        None, "--alignments", "-a", dir_okay=False, exists=True, help="Path to sequence alignments in fasta format"
    ),
    tree: Optional[Path] = typer.Option(None, "--tree", "-t", dir_okay=False, exists=True, help="Path to tree file"),
):
    args = [testfile]
    if alignments is not None:
        args.extend(["--alignments", alignments])
    if tree is not None:
        args.extend(["--tree", tree])

    pytest.main(args)


app = typer.run(cli)
