import copy
import re
from datetime import datetime
from io import StringIO
from typing import Dict, List, Optional, Union
from warnings import warn

from Bio import Phylo as Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo.BaseTree import Tree as BioTree
from dateutil.parser import parse
from treetime import GTR, TreeTime
from treetime.utils import DateConversion, datetime_from_numeric, numeric_date

from ..utils import PhytestObject, PhytestWarning, assert_or_warn, default_date_patterns


class Tree(PhytestObject, BioTree):
    @classmethod
    def read(cls, tree_path, tree_format) -> 'Tree':
        tree = Phylo.read(tree_path, tree_format)
        return cls(root=tree.root, rooted=tree.rooted, id=tree.id, name=tree.name)

    @classmethod
    def parse(cls, tree_path, tree_format) -> 'Tree':
        trees = Phylo.parse(tree_path, tree_format)
        return (cls(root=tree.root, rooted=tree.rooted, id=tree.id, name=tree.name) for tree in trees)

    @classmethod
    def read_str(cls, tree_str: str, tree_format: str = "newick") -> 'Tree':
        data = StringIO(tree_str)
        return cls.read(data, tree_format)

    def parse_tip_dates(
        self,
        *,
        patterns=None,
        date_format: Optional[str] = None,
        decimal_year: bool = False,
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
        warning: bool = False,
    ):
        """
        Asserts that that the number of tips meets the specified criteria.

        Args:
            tips (int, optional): If set, then number of tips must be equal to this value. Defaults to None.
            min (int, optional): If set, then number of tips must be equal to or greater than this value. Defaults to None.
            max (int, optional): If set, then number of tips must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        number_of_tips = len(self.get_terminals())
        if tips is not None:
            assert_or_warn(
                number_of_tips == tips,
                warning,
                f"The number of tips ({number_of_tips}) which is different from the required number of tips ({tips}).",
            )
        if min is not None:
            assert_or_warn(
                number_of_tips >= min,
                warning,
                f"The number of tips ({number_of_tips}) is less than the minimum ({min}).",
            )
        if max is not None:
            assert_or_warn(
                number_of_tips <= max,
                warning,
                f"The number of tips ({number_of_tips}) is greater than the maximum ({max}).",
            )

    def assert_is_bifurcating(self, *, warning: bool = False):
        """
        Asserts that the tree is bifurcating.

        The root may have 3 descendents and still be considered part of a bifurcating tree, because it has no ancestor.

        Args:
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        assert_or_warn(
            self.is_bifurcating(),
            warning,
            f"The tree is not bifurcating.",
        )

    def assert_is_monophyletic(self, tips: List[Clade], *, warning: bool = False):
        """
        Asserts that the specified tips form a monophyletic group.

        Args:
            tips (List[Clade]): List of terminal nodes (tips).
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        assert_or_warn(
            self.is_monophyletic(tips),
            warning,
            f"The group \'{', '.join([tip.name for tip in tips])}\' is paraphyletic!",
        )

    def assert_total_branch_length(
        self,
        length: Optional[float] = None,
        *,
        min: Optional[float] = None,
        max: Optional[float] = None,
        warning: bool = False,
    ):
        """
        Asserts that that the total brach length meets the specified criteria.

        Args:
            length (float, optional): If set, then total brach length must be equal to this value. Defaults to None.
            min (float, optional): If set, then total brach length must be equal to or greater than this value. Defaults to None.
            max (float, optional): If set, then total brach length must be equal to or less than this value. Defaults to None.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
        """
        total_branch_length = self.total_branch_length()
        if length is not None:
            assert_or_warn(
                total_branch_length == length,
                warning,
                f"The total branch length ({total_branch_length}) is not equal to the required length ({length}).",
            )
        if min is not None:
            assert_or_warn(
                total_branch_length >= min,
                warning,
                f"The total branch length ({total_branch_length}) is less than the minimum ({length}).",
            )
        if max is not None:
            assert_or_warn(
                total_branch_length <= max,
                warning,
                f"The total branch length ({total_branch_length}) is greater than the maximum ({length}).",
            )

    def assert_tip_regex(
        self,
        patterns: Union[List[str], str],
        *,
        warning: bool = False,
    ):
        """
        Asserts that all the tips match at least one of a list of regular expression patterns.

        Args:
            patterns (Union[List[str], str]): The regex pattern(s) to match to.
                If a string, then every tip must match that pattern.
                If a list then each tip must match at least one pattern in the list.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
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
            assert_or_warn(
                matches,
                warning,
                f"Tip {tip.name} does not match any of the regex patterns in: '{patterns}'.",
            )

    def root_to_tip_regression(
        self,
        *,
        dates: Optional[Dict] = None,
        alignment: Optional[MultipleSeqAlignment] = None,
        sequence_length: Optional[int] = None,
        clock_filter: float = 3.0,
        gtr: Union[GTR, str] = 'JC69',
        root_method: str = 'least-squares',
        allow_negative_rate: bool = False,
        keep_root: bool = False,
        covariation: bool = False,
    ):
        """
        Performs a root-to-tip regression to determine how clock-like a tree is.

        Args:
            dates (Dict, optional): The tip dates as a dictionary with the tip name as the key and the date as the value.
                If not set, then it parses the tip dates to generate this dictionary using the `parse_tip_dates` method.
            alignment (MultipleSeqAlignment, optional): The alignment associated with this tree. Defaults to None.
            sequence_length (int, optional): The sequence length of the alignment. Defaults to None.
            clock_filter (float, optional): The number of interquartile ranges from regression beyond which to ignore.
                This provides a way to ignore tips that don't follow a loose clock.
                Defaults to 3.0.
            gtr (GTR, str, optional): The molecular evolution model. Defaults to 'JC69'.
            allow_negative_rate (bool, optional): Whether or not a negative clock rate is allowed.
                For trees with little temporal signal, it can be set to True to achieve essentially mid-point rooting.
                Defaults to False.
            keep_root (bool, optional): Keeps the current root of the tree.
                If False, then a new optimal root is sought. Defaults to False.
            root_method (str, optional): The method used to reroot the tree if `keep_root` is False.
                Valid choices are: 'min_dev', 'least-squares', and 'oldest'.
                Defaults to 'least-squares'.
            covariation (bool, optional): Accounts for covariation when estimating rates or rerooting. Defaults to False.
        """
        dates = dates or self.parse_tip_dates()

        # Convert datetimes to floats with decimal years if necessary
        dates = {name: numeric_date(date) if isinstance(date, datetime) else date for name, date in dates.items()}

        regression = TreeTime(
            dates=dates,
            tree=copy.deepcopy(self),
            aln=alignment,
            gtr=gtr,
            seq_len=sequence_length,
        )

        if clock_filter:
            bad_nodes = [node.name for node in regression.tree.get_terminals() if node.bad_branch]
            regression.clock_filter(n_iqd=clock_filter, reroot=root_method or 'least-squares')
            bad_nodes_after = [node.name for node in regression.tree.get_terminals() if node.bad_branch]
            if len(bad_nodes_after) > len(bad_nodes):
                warn(
                    "The following leaves don't follow a loose clock and "
                    "will be ignored in rate estimation:\n\t" + "\n\t".join(set(bad_nodes_after).difference(bad_nodes)),
                    PhytestWarning,
                )

        if not keep_root:
            if covariation:  # this requires branch length estimates
                regression.run(root="least-squares", max_iter=0, use_covariation=covariation)

            assert root_method in ['min_dev', 'least-squares', 'oldest']
            regression.reroot(root_method, force_positive=not allow_negative_rate)

        regression.get_clock_model(covariation=covariation)
        return regression

    def assert_root_to_tip(
        self,
        *,
        regression: Optional[TreeTime] = None,
        min_r_squared: Optional[float] = None,
        min_rate: Optional[float] = None,
        max_rate: Optional[float] = None,
        min_root_date: Optional[float] = None,
        max_root_date: Optional[float] = None,
        assert_valid_confidence: bool = False,
        warning: bool = False,
        **kwargs,
    ):
        """
        Performs a root-to-tip regression to determine how clock-like a tree is.

        Args:
            regression (TreeTime, optional): The root-to-tip regression for this tree.
                If None, then this regression is calculated using the `root_to_tip_regression` method.
            min_r_squared (float, optional): If set, then R^2 must be equal or greater than this value. Defaults to None.
            min_rate (float, optional): If set, then the clock rate must be equal or greater than this value. Defaults to None.
            max_rate (float, optional): If set, then the clock rate must be equal or less than this value. Defaults to None.
            min_root_date (float, optional): If set, then the interpolated root date must be equal or greater than this value. Defaults to None.
            max_root_date (float, optional): If set, then the interpolated root date must be equal or less than this value. Defaults to None.
            assert_valid_confidence (bool, optional): If set then `valid_confidence` in the regression must be True. Defaults to False.
            warning (bool): If True, raise a warning instead of an exception. Defaults to False.
                This flag can be set by running this method with the prefix `warn_` instead of `assert_`.
            **kwargs: Keyword arguments for the `root_to_tip_regression` method.
        """
        regression = regression or self.root_to_tip_regression(**kwargs)
        clock_model = DateConversion.from_regression(regression.clock_model)
        root_date = clock_model.numdate_from_dist2root(0.0)

        if min_r_squared is not None:
            assert_or_warn(
                clock_model.r_val**2 >= min_r_squared,
                warning,
                f"The R-squared value from the root-to-tip regression '{clock_model.r_val**2}' "
                "is less than the minimum allowed R-squarred '{min_r_squared}'.",
            )

        if min_rate is not None:
            assert_or_warn(
                clock_model.clock_rate >= min_rate,
                warning,
                f"Inferred clock rate '{clock_model.clock_rate}' is less than the minimum allowed clock rate '{min_rate}'.",
            )

        if max_rate is not None:
            assert_or_warn(
                clock_model.clock_rate <= max_rate,
                warning,
                f"Inferred clock rate '{clock_model.clock_rate}' is greater than the maximum allowed clock rate '{max_rate}'.",
            )

        if min_root_date is not None:
            assert_or_warn(
                root_date >= min_root_date,
                warning,
                f"Inferred root date '{root_date}' is less than the minimum allowed root date '{min_root_date}'.",
            )

        if max_root_date is not None:
            assert_or_warn(
                root_date <= max_root_date,
                warning,
                f"Inferred root date '{root_date}' is greater than the maximum allowed root date: '{max_root_date}'.",
            )

        if assert_valid_confidence:
            assert_or_warn(
                clock_model.valid_confidence, warning, f"The `clock_model.valid_confidence` variable is False."
            )
