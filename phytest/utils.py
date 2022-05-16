from warnings import warn


class PhytestWarning(Warning):
    pass


class PhytestAssertion(AssertionError):
    pass


def assert_or_warn(statement, warning, message):
    if statement:
        return

    if warning:
        warn(message, PhytestWarning)
    else:
        raise PhytestAssertion(message)


def default_date_patterns():
    return [
        r"\d{4}\.?\d*$",
        r"\d{4}-\d{2}-\d{2}",
    ]
