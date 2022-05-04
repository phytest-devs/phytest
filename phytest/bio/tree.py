import re
import copy
from dateutil.parser import parse
from datetime import datetime
from io import StringIO
from typing import Union, Optional, Dict, List
from warnings import warn

from Bio import Phylo as Phylo
from Bio.Phylo.BaseTree import Tree as BioTree
from Bio.Align import MultipleSeqAlignment

from treetime import TreeTime, GTR
from treetime.utils import numeric_date, datetime_from_numeric, DateConversion
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

    def assert_root_to_tip(
        self,
        *,
        dates: Optional[Dict] = None,
        alignment: Optional[MultipleSeqAlignment] = None,
        sequence_length: Optional[int] = None,
        clock_filter:float = 3.0,
        gtr:Union[GTR,str] = 'JC69',
        root_method:str = 'least-squares',
        allow_negative_rate:bool = False,
        keep_root:bool = False,
        covariation:bool = False,
        min_r_squared: Optional[float] = None,
        min_rate: Optional[float] = None,
        max_rate: Optional[float] = None,
        assert_valid_confidence: bool = False,
    ):
        """
        Performs a root-to-tip regression to determine how clock-like a tree is.

        Args:
            dates (Optional[Dict], optional): The tip dates as a dictionary with the tip name as the key and the date as the value. 
                If not set, then it parses the tip dates to generate this dictionary using the `parse_tip_dates` method.
            alignment (Optional[MultipleSeqAlignment], optional): The alignment associated with this tree. Defaults to None.
            sequence_length (Optional[int], optional): The sequence length of the alignment. Defaults to None.
            clock_filter (float, optional): The number of interquartile ranges from regression beyond which to ignore.
                This provides a way to ignore tips that don't follow a loose clock. 
                Defaults to 3.0.
            gtr (GTR, str, optional): The molecular evolution model. Defaults to 'JC69'.
            allow_negative_rate (bool, optional): Whether or not a negative clock rate is allowed. 
                For trees with little temporal signal, it can be set to True to achieve essentially mid-point rooting.
                Defaults to False.
            keep_root (bool, optional): Keeps the current root of the tree. If False, then a new optimal root is . Defaults to False.
            root_method (str, optional): The method used to reroot the tree if `keep_root` is False. 
                Valid choices are: 'min_dev', 'least-squares', and 'oldest'.
                Defaults to 'least-squares'.
            covariation (bool, optional): Account for covariation when estimating rates or rerooting. Defaults to False.
            min_r_squared (Optional[float], optional): If set, then R^2 must be equal or greater than this value. Defaults to None.
            min_rate (Optional[float], optional): If set, then the clock rate must be equal or greater than this value. Defaults to None.
            max_rate (Optional[float], optional): If set, then the clock rate must be equal or lesser than this value. Defaults to None.
            assert_valid_confidence (bool, optional): If set then `valid_confidence` in the regression must be True. Defaults to False.
        """
        dates = dates or self.parse_tip_dates()

        # Convert datetimes to floats with decimal years if necessary
        dates = {name:numeric_date(date) if isinstance(date, datetime) else date for name,date in dates.items()}

        treetime = TreeTime(
            dates=dates, 
            tree=copy.deepcopy(self), 
            aln=alignment, 
            gtr=gtr,
            seq_len=sequence_length,
        )

        if clock_filter:
            bad_nodes = [node.name for node in treetime.tree.get_terminals() if node.bad_branch]
            treetime.clock_filter(n_iqd=clock_filter, reroot=root_method or 'least-squares')
            bad_nodes_after = [node.name for node in treetime.tree.get_terminals() if node.bad_branch]
            if len(bad_nodes_after) > len(bad_nodes):
                print(
                    "The following leaves don't follow a loose clock and "
                    "will be ignored in rate estimation:\n\t"
                    +"\n\t".join(set(bad_nodes_after).difference(bad_nodes))
                )

        if not keep_root:
            if covariation: # this requires branch length estimates
                treetime.run(root="least-squares", max_iter=0, use_covariation=covariation)

            assert root_method in ['min_dev', 'least-squares', 'oldest']
            treetime.reroot(root_method, force_positive=not allow_negative_rate)

        treetime.get_clock_model(covariation=covariation)
        clock_model = DateConversion.from_regression(treetime.clock_model)

        if min_r_squared is not None:
            assert clock_model.r_val**2 >= min_r_squared

        if min_rate is not None:
            assert clock_model.clock_rate >= min_rate

        if max_rate is not None:
            assert clock_model.clock_rate <= max_rate

        if assert_valid_confidence:
            assert clock_model.valid_confidence
