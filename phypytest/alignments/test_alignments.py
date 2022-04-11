import re

regex_invalid = re.compile(r"[^AGTCN-]")


def test_invalid_base(sequence, should_fix):
    seq_id, seq = sequence
    assert not regex_invalid.search(seq), f"Invalid base found in \'{seq_id.lstrip('>')}\'."
    if should_fix:
        print("No worries chief... I can fix that! 🛠 ")
