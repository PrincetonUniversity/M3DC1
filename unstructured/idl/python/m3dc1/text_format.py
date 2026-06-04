from __future__ import annotations

import re


def idl_to_plain(text: str) -> str:
    """
    Convert common IDL plot control tokens to plain text.
    This avoids shell/matplotlib issues with strings like '!X', '!6', '!D', '!N'.
    """
    s = str(text)
    # Drop IDL font/format codes such as !6, !8, !9, !X, !N, !D, !U, !3, etc.
    s = re.sub(r"![A-Za-z0-9]+", "", s)
    # Collapse repeated whitespace left behind by removed tokens.
    s = re.sub(r"\s+", " ", s).strip()
    return s
