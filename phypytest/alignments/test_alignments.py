import re

regex_invalid = re.compile(r"[^AGTCN-]")


def test_invalid_base(sequence):
    seq_id, seq = sequence
    assert not regex_invalid.search(seq), f"Invalid base found in \'{seq_id.lstrip('>')}\'."
