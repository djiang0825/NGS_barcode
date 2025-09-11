"""
Unified CLI for Barcadia.
Usage:
  barcadia generate [options...]   -> delegates to barcadia.generate_barcodes.main(argv)
  barcadia validate [options...]   -> delegates to barcadia.validate_barcodes.main(argv)
"""

import sys
from . import generate_barcodes as gen
from . import validate_barcodes as val

TOP_USAGE = (
    "Usage:\n"
    "  barcadia generate [options...]\n"
    "  barcadia validate [options...]\n"
    "\n"
    "Use `barcadia <subcommand> --help` to see subcommand options.\n"
)

def main() -> int:
    # No subcommand → show top-level help
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help"}:
        print(TOP_USAGE, file=sys.stderr)
        return 0

    cmd, argv = sys.argv[1], sys.argv[2:]

    if cmd == "generate":
        # gen.main must accept argv: list[str] | None
        return gen.main(argv) or 0

    if cmd == "validate":
        # val.main must accept argv: list[str] | None
        return val.main(argv) or 0

    # Unknown subcommand
    print(f"Unknown subcommand: {cmd}\n\n{TOP_USAGE}", file=sys.stderr)
    return 2

if __name__ == "__main__":
    sys.exit(main())
