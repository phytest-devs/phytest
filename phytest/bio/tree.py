import re
from dateutil.parser import parse
from io import StringIO
from typing import Union, Optional, Dict, List
from warnings import warn

from Bio import Phylo as Phylo
from Bio.Phylo.BaseTree import Tree as BioTree
from Bio.Align import MultipleSeqAlignment

from treetime.utils import numeric_date, datetime_from_numeric
from ..utils import default_date_patterns

class Tree(BioTree):
    @classmethod
    def read(cls, tree_path, tree_format) -> 'Tree':
        tree = Phylo.read(tree_path, tree_format)
        return Tree(root=tree.root, rooted=tree.rooted, id=tree.id, name=tree.name)

    @classmethod
    def read_str(cls, tree_str:str, tree_format:str="newick") -> 'Tree':
        data = StringIO(tree_str)
        return cls.read( data, tree_format )

    def parse_tip_dates(
        self, 
        *,
        patterns = None,        
        date_format: Optional[str] = None,
        decimal_year:bool = False,
    ):
        patterns = patterns or default_date_patterns()
        if isinstance(patterns, str):
            patterns = [patterns]

        dates = {}

        compiled_patterns = [re.compile(pattern_string) for pattern_string in patterns]
        for tip in self.find_elements(terminal=True):
            for pattern in compiled_patterns:
                m = pattern.search(tip.name)
                if m:
                    matched_str = m.group(0)
                    if re.match(r"^\d+\.?\d*$", matched_str):
                        date = datetime_from_numeric(float(matched_str))
                    else:
                        date = parse(matched_str, date_format)
                    
                    dates[tip.name] = date
                    break
        
        if decimal_year:
            dates = {key: numeric_date(value) for key, value in dates.items()}
        
        return dates

    def assert_number_of_tips(
        self,
        tips: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ):
        number_of_tips = len(self.get_terminals())
        if tips is not None:
            assert number_of_tips == tips
        if min is not None:
            assert number_of_tips >= min
        if max is not None:
            assert number_of_tips <= max
        if warning is not None and number_of_tips != warning:
            warn(f"Number of tips '{number_of_tips}' != {warning}")

    def assert_is_bifurcating(self):
        assert self.is_bifurcating()

    def assert_total_branch_length(
        self,
        length: Optional[int] = None,
        *,
        min: Optional[int] = None,
        max: Optional[int] = None,
        warning: Optional[int] = None,
    ):
        total_branch_length = self.total_branch_length()
        if length is not None:
            assert total_branch_length == length
        if min is not None:
            assert total_branch_length >= min
        if max is not None:
            assert total_branch_length <= max
        if warning is not None and total_branch_length != warning:
            warn(f"Total branch length '{total_branch_length}' != {warning}")

    def assert_tip_regex(
        self,
        patterns: Union[List[str], str],
    ):
        """
        Asserts that all the tips match at least one of a list of regular expression patterns.

        Args:
            patterns (Union[List[str], str]): The regex pattern(s) to match to. 
                If a string, then every tip must match that pattern. 
                If a list then each tip must match at least one pattern in the list.
        """
        if isinstance(patterns, str):
            patterns = [patterns]

        compiled_patterns = [re.compile(pattern_string) for pattern_string in patterns]

        for tip in self.find_elements(terminal=True):
            matches = False
            for pattern in compiled_patterns:
                if pattern.search(tip.name):
                    matches = True
                    break
            assert matches, f"Tip {tip.name} does not match any of the regex patterns in: '{patterns}'."

