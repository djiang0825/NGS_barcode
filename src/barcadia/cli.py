"""
Unified CLI for Barcadia.
Usage:
  barcadia generate [options...]   -> delegates to barcadia.generate_barcodes.main(argv)
  barcadia validate [options...]   -> delegates to barcadia.validate_barcodes.main(argv)
"""

import sys
from importlib.metadata import version

from . import generate_barcodes as gen
from . import validate_barcodes as val

TOP_USAGE = (
    "Barcadia - A high-performance, memory-efficient toolkit for fast generation and validation of large-scale NGS barcodes\n"
    "\n"
    "Usage:\n"
    "  barcadia <command> [options...]\n"
    "\n"
    "Commands:\n"
    "  generate    Generate high-performance DNA barcodes for NGS applications\n"
    "  validate    Validate DNA barcodes against quality filters\n"
    "\n"
    "Examples:\n"
    "  barcadia --help\n"
    "  barcadia generate --help\n"
    "  barcadia validate --help\n"
    "  barcadia generate --count 1000 --length 12\n"
    "  barcadia validate --input test/barcodes.txt\n"
    "\n"
    "Global options:\n"
    "  --help, -h     Show this help message\n"
    "  --version, -v  Show version information\n"
)


def main() -> int:
    # Handle version flag
    if len(sys.argv) >= 2 and sys.argv[1] in {"-v", "--version"}:
        print(version("barcadia"))
        return 0

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
