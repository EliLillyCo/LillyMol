#!/usr/bin/env python3

# Check the status of all the .qry files in the directory and see whether each query is in
# one of the --list-file or not.

# python check_query_file_status.py --query-dir . --list-file demerits --list-file reject1 --list-file reject2

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Check that every .qry file in a directory appears in exactly one "
            "query-list file."
        )
    )
    parser.add_argument(
        "--query-dir",
        type=Path,
        required=True,
        help="Directory containing .qry files.",
    )
    parser.add_argument(
        "--list-file",
        type=Path,
        action="append",
        required=True,
        dest="list_files",
        help=(
            "File containing query file names. Specify once per list file. "
            "Expected three times, but any number is accepted."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional report file. If omitted, report is written to stdout.",
    )
    return parser.parse_args()


def query_files_in_directory(query_dir: Path) -> set[str]:
    if not query_dir.is_dir():
        raise ValueError(f"Not a directory: {query_dir}")

    return {path.name for path in query_dir.glob("*.qry") if path.is_file()}


def read_list_file(path: Path) -> list[str]:
    if not path.is_file():
        raise ValueError(f"List file does not exist: {path}")

    result: list[str] = []

    with path.open("r", encoding="utf-8") as reader:
        for line_number, line in enumerate(reader, start=1):
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            # Allow comments after entries:
            #   foo.qry   # comment
            if "#" in line:
                line = line.split("#", 1)[0].strip()

            if not line:
                continue

            # Normalize ./foo.qry or /some/path/foo.qry to foo.qry.
            name = Path(line).name

            if not name.endswith(".qry"):
                raise ValueError(
                    f"{path}:{line_number}: entry does not end in .qry: {line!r}"
                )

            result.append(name)

    return result


def build_membership(list_files: list[Path]) -> dict[str, list[str]]:
    membership: dict[str, list[str]] = defaultdict(list)

    for list_file in list_files:
        list_name = list_file.name
        entries = read_list_file(list_file)

        for entry in entries:
            membership[entry].append(list_name)

    return membership


def make_report(
    *,
    query_files: set[str],
    membership: dict[str, list[str]],
) -> str:
    lines: list[str] = []

    lines.append("Query file mapping")
    lines.append("==================")
    lines.append("")

    for query in sorted(query_files):
        owners = membership.get(query, [])

        if len(owners) == 0:
            owner_text = "MISSING"
        elif len(owners) == 1:
            owner_text = owners[0]
        else:
            owner_text = "DUPLICATE: " + ", ".join(owners)

        lines.append(f"{query}\t{owner_text}")

    lines.append("")
    lines.append("Missing from all list files")
    lines.append("===========================")
    lines.append("")

    missing = sorted(query for query in query_files if len(membership.get(query, [])) == 0)

    if missing:
        lines.extend(missing)
    else:
        lines.append("None")

    lines.append("")
    lines.append("Present in more than one list file")
    lines.append("==================================")
    lines.append("")

    duplicates = sorted(
        query for query in query_files if len(membership.get(query, [])) > 1
    )

    if duplicates:
        for query in duplicates:
            lines.append(f"{query}\t{', '.join(membership[query])}")
    else:
        lines.append("None")

    lines.append("")
    lines.append("Listed but not found in query directory")
    lines.append("======================================")
    lines.append("")

    listed_queries = set(membership)
    extra_listed = sorted(listed_queries - query_files)

    if extra_listed:
        for query in extra_listed:
            lines.append(f"{query}\t{', '.join(membership[query])}")
    else:
        lines.append("None")

    lines.append("")
    lines.append("Summary")
    lines.append("=======")
    lines.append("")
    lines.append(f"Total .qry files in directory: {len(query_files)}")
    lines.append(f"Missing from all list files:   {len(missing)}")
    lines.append(f"Duplicated across list files:  {len(duplicates)}")
    lines.append(f"Listed but not found:          {len(extra_listed)}")

    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()

    query_files = query_files_in_directory(args.query_dir)
    membership = build_membership(args.list_files)

    report = make_report(
        query_files=query_files,
        membership=membership,
    )

    if args.output is None:
        print(report, end="")
    else:
        args.output.write_text(report, encoding="utf-8")

    missing_count = sum(1 for query in query_files if len(membership.get(query, [])) == 0)
    duplicate_count = sum(1 for query in query_files if len(membership.get(query, [])) > 1)
    listed_but_not_found_count = len(set(membership) - query_files)

    if missing_count or duplicate_count or listed_but_not_found_count:
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
