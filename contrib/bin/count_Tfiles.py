#!/usr/bin/env python3

import argparse
import logging
from collections import defaultdict
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sum counts for tokens across two-column files."
    )
    parser.add_argument(
        "files",
        nargs="+",
        type=Path,
        help="Input files. Each line should contain: count token",
    )
    parser.add_argument(
        "--integer",
        action="store_true",
        help="Parse and print counts as integers instead of floats.",
    )

    args = parser.parse_args()

    totals = defaultdict(float)

    for path in args.files:
        with path.open("r", encoding="utf-8") as reader:
            for line_number, line in enumerate(reader, start=1):
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                fields = line.split(maxsplit=1)

                if len(fields) != 2:
                    logging.warning("Skipping malformed line %s", line)
                    continue
                    # raise ValueError(
                    #     f"{path}:{line_number}: expected 2 columns, got {len(fields)}"
                    # )

                count_text, token = fields

                try:
                    count = int(count_text) if args.integer else float(count_text)
                except ValueError as err:
                    raise ValueError(
                        f"{path}:{line_number}: invalid count: {count_text}"
                    ) from err

                totals[token] += count

    for token in sorted(totals):
        total = totals[token]

        if args.integer:
            total = int(total)

        print(token, total)


if __name__ == "__main__":
    main()
