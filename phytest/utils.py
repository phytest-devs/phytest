from functools import partial
from typing import List
from warnings import warn


class PhytestWarning(Warning):
    pass


class PhytestAssertion(AssertionError):
    pass


def assert_or_warn(statement, warning, *messages):
    if statement:
        return

    message = "\n".join(messages)
    if warning:
        warn(message, PhytestWarning)
    else:
        raise PhytestAssertion(message)


def default_date_patterns():
    return [
        r"\d{4}\.?\d*$",
        r"\d{4}-\d{2}-\d{2}",
    ]


class PhytestObject:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Add partial methods with the warning flag set to True
        for method_name in self.assertion_method_names():
            method = getattr(self, method_name)
            truncated_name = method_name[len("assert") :]
            warning_name = f"warn{truncated_name}"
            setattr(self, warning_name, partial(method, warning=True))

    def assertion_method_names(self) -> List[str]:
        """
        Returns a list with the names of the methods used to make assertion statements.
        """
        return [
            attribute
            for attribute in dir(self)
            if callable(getattr(self, attribute)) and attribute.startswith("assert_")
        ]
