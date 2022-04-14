import inspect
import os
from pathlib import Path
from typing import Optional
from xmlrpc.client import Boolean

import pytest


def main(
    testfile: Optional[Path] = None,
    alignment: Optional[Path] = None,
    alignment_format: Optional[str] = 'fasta',
    tree: Optional[Path] = None,
    tree_format: Optional[str] = 'newick',
    verbose: bool = False,
    report: bool = False,
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
    if report:
        args.extend(["--html=report.html", "--self-contained-html"])
    pytest.main(args, plugins=['phytest'])
