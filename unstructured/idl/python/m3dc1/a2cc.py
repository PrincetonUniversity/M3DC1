from __future__ import annotations

from pathlib import Path
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .readaeqdsk import AeqdskData


def a2cc(filename: str | Path) -> AeqdskData:
    # program a2cc
    #
    #   call load_eqdsk_a(filename)
    #
    # Python port:
    #   return the parsed a-file data structure produced by the readaeqdsk module.
    from m3dc1 import load_eqdsk_a

    return load_eqdsk_a(filename)


def _main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("Usage: a2cc <aeqdsk>", file=sys.stderr)
        return 1
    a2cc(args[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
