import re
from warnings import warn

regex_invalid = re.compile(r"[^AGTCN-]")


def test_invalid_base(sequence, should_fix):
    seq_id, seq = sequence
    assert not regex_invalid.search(seq), f"Invalid base found in \'{seq_id.lstrip('>')}\'."
    if should_fix:
        print("No worries chief... I can fix that! 🛠 ")


def test_width(sequence, width=100):
    assert len(sequence[1].strip()) == width


def test_count_Ns(sequence, max=1, warning=0):
    ns = sequence.count("N")
    if ns > warning:
        warn(f"Found {ns} N bases (> {warning})")
    assert ns < max
