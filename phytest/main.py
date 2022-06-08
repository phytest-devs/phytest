import inspect
import os
from pathlib import Path
from typing import Optional

import pytest


def main(
    testfile: Optional[Path] = None,
    alignment: Optional[Path] = None,
    alignment_format: Optional[str] = 'fasta',
    tree: Optional[Path] = None,
    tree_format: Optional[str] = 'newick',
    data: Optional[Path] = None,
    verbose: bool = False,
    report: Optional[Path] = None,
    expression: str = None,
):
    if not testfile:
        testfile = Path(os.path.abspath((inspect.stack()[1])[1]))
    args = [testfile]
    if not verbose:
        args.extend(["-ra", "--tb=no", "--no-header"])
    if alignment is not None:
        args.extend(["--alignment", alignment])
        args.extend(["--alignment-format", alignment_format])
    if tree is not None:
        args.extend(["--tree", tree])
        args.extend(["--tree-format", tree_format])
    if data is not None:
        args.extend(["--data", data])
    if report:
        if not str(report).endswith('.html'):
            raise ValueError(f"Report must use .html extension.")
        args.extend([f"--html={report}", "--self-contained-html"])
    if expression:
        args.extend(["-k", expression])
    exit_code = pytest.main(args, plugins=['phytest'])
    return exit_code
