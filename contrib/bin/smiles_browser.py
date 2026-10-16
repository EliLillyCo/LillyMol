#!/usr/bin/env python3

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Draw

import argparse
import csv
from contextlib import nullcontext
import html
import json
import re
import sys
from dataclasses import dataclass
from typing import Optional
from urllib.parse import quote_plus


@dataclass
class LinkSpec:
    column: int
    template: str


@dataclass
class ColorRule:
    op: str
    value: str
    color: str
    regex: Optional[re.Pattern] = None


def mol_to_svg(smiles: str, sanitize=False, width=250, height=180) -> str:
    mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
    if mol is None:
        return f"<span class='bad'>{html.escape(smiles)} bad</span>"

    return Draw.MolsToGridImage(
        [mol],
        molsPerRow=1,
        subImgSize=(width, height),
        useSVG=True,
    )


def separator_from_args(args):
    if args.separator:
        if args.separator == r"\t":
            return "\t"
        return args.separator

    if args.input.lower().endswith(".csv"):
        return ","

    return None   # whitespace split


def parse_smiles_columns(values) -> list[int]:
    columns = []

    def flattened(items):
        for item in items or ["1"]:
            if isinstance(item, (list, tuple)):
                yield from flattened(item)
            else:
                yield item

    for value in flattened(values):
        for token in str(value).split(","):
            token = token.strip()
            if not token:
                continue
            try:
                column = int(token)
            except ValueError:
                sys.exit(f"--smiles_column: invalid column '{token}'")
            if column < 1:
                sys.exit("--smiles_column values must be positive 1-based column numbers")
            columns.append(column - 1)

    if not columns:
        sys.exit("--smiles_column must specify at least one column")

    seen = set()
    result = []
    for column in columns:
        if column not in seen:
            seen.add(column)
            result.append(column)

    return result


def split_smiles_and_extra(tokens: list[str], smiles_columns: list[int], record_number: int):
    fields = [value.strip() for value in tokens]
    smiles = []
    for column in smiles_columns:
        if column >= len(fields):
            sys.exit(f"Input record {record_number} does not have SMILES column {column + 1}")
        smiles.append(fields[column])

    smiles_set = set(smiles_columns)
    extra = [value for i, value in enumerate(fields) if i not in smiles_set]
    return smiles, extra, fields


def read_smiles_file(fname, separator=None, has_header=False, max_molecules=None, smiles_columns=None):
    rows = []
    max_extra = 0
    headers = None
    display_headers = None
    max_fields = 0
    smiles_columns = smiles_columns or [0]

    # Leave stdin open; only streams opened here belong to this reader.
    source = nullcontext(sys.stdin) if fname == "-" else open(fname, newline="")
    with source as f:
        if separator == ",":
            reader = csv.reader(f)
            has_header = True
        elif separator is not None:
            reader = csv.reader(f, delimiter=separator)
        else:
            reader = None

        if reader is not None:
            for record_number, tokens in enumerate(reader, 1):
                if not tokens:
                    continue

                if has_header and headers is None:
                    _, headers, display_headers = split_smiles_and_extra(tokens, smiles_columns, record_number)
                    continue

                smiles, extra, fields = split_smiles_and_extra(tokens, smiles_columns, record_number)

                if not any(smiles):
                    continue

                rows.append((record_number, smiles, extra, fields))
                max_extra = max(max_extra, len(extra))
                max_fields = max(max_fields, len(fields))

                if max_molecules is not None and len(rows) >= max_molecules:
                    break
        else:
            for record_number, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                tokens = line.split()

                if has_header and headers is None:
                    _, headers, display_headers = split_smiles_and_extra(tokens, smiles_columns, record_number)
                    continue

                smiles, extra, fields = split_smiles_and_extra(tokens, smiles_columns, record_number)

                rows.append((record_number, smiles, extra, fields))
                max_extra = max(max_extra, len(extra))
                max_fields = max(max_fields, len(fields))

                if max_molecules is not None and len(rows) >= max_molecules:
                    break

    smiles_set = set(smiles_columns)
    if display_headers is None:
        display_headers = []
        field_number = 1
        structure_number = 1
        for column in range(max_fields):
            if column in smiles_set:
                if len(smiles_columns) == 1:
                    display_headers.append("Structure")
                else:
                    display_headers.append(f"Structure {structure_number}")
                structure_number += 1
            else:
                display_headers.append(f"Field {field_number}")
                field_number += 1

    if headers is None:
        headers = [display_headers[i] for i in range(len(display_headers)) if i not in smiles_set]
    elif len(headers) < max_extra:
        headers.extend(f"Field {i + 1}" for i in range(len(headers), max_extra))

    return rows, display_headers, headers


def column_from_string(column: str, headers: list[str], option_name: str) -> int:
    """Return a zero-based extra-field column index from a 1-based number or header."""
    try:
        value = int(column)
    except ValueError:
        try:
            return headers.index(column)
        except ValueError:
            sys.exit(f"{option_name}: no column named '{column}'")

    if value < 1 or value > len(headers):
        sys.exit(f"{option_name}: column {value} is outside the valid range 1..{len(headers)}")

    return value - 1


def parse_link_columns(link_columns: list[str], headers: list[str]) -> list[LinkSpec]:
    result = []

    for item in link_columns or []:
        if "=" not in item:
            sys.exit("--link-column must be of the form column=url_template")

        column, template = item.split("=", 1)
        column = column.strip()
        template = template.strip()

        if not column or not template:
            sys.exit("--link-column must have both a column and a url_template")

        result.append(LinkSpec(column_from_string(column, headers, "--link-column"), template))

    return result


def parse_image_columns(image_columns: list[str], headers: list[str]) -> set[int]:
    return {column_from_string(column.strip(), headers, "--image_column")
            for column in image_columns or []}


def parse_color_rules(color_rules: list[str]) -> list[ColorRule]:
    result = []

    for item in color_rules or []:
        parts = item.split(":", 2)
        if len(parts) != 3:
            sys.exit("--color-rule must be of the form op:value:css_color")

        op, value, color = [x.strip() for x in parts]
        if op not in {"lt", "le", "eq", "ne", "ge", "gt", "contains", "startswith", "regex"}:
            sys.exit(f"--color-rule: unrecognised operator '{op}'")

        compiled = None
        if op == "regex":
            try:
                compiled = re.compile(value)
            except re.error as e:
                sys.exit(f"--color-rule: invalid regex '{value}': {e}")

        result.append(ColorRule(op, value, color, compiled))

    return result


def link_for_value(template: str, value: str) -> str:
    """Build a URL from a template.

    The template can contain {}, {value}, {raw}, or {quoted}. The {} and
    {quoted} forms URL-encode the value. {value} and {raw} insert it unchanged.
    """
    quoted = quote_plus(value)

    if "{}" in template:
        return template.replace("{}", quoted)

    try:
        return template.format(value=value, raw=value, quoted=quoted)
    except (IndexError, KeyError, ValueError):
        # If the template contains braces that are not meant for Python format,
        # fall back to appending the quoted value.
        return template + quoted


def maybe_link_value(value: str, column: int, link_specs: list[LinkSpec]) -> str:
    escaped_value = html.escape(value)
    if not value:
        return escaped_value

    for spec in link_specs:
        if spec.column == column:
            url = html.escape(link_for_value(spec.template, value), quote=True)
            return f'<a href="{url}" target="_blank" rel="noopener noreferrer">{escaped_value}</a>'

    return escaped_value


def image_cell(value: str, image_width: int) -> str:
    if not value:
        return ""

    escaped_value = html.escape(value)
    escaped_src = html.escape(value, quote=True)
    return (
        f'<a href="{escaped_src}" target="_blank" rel="noopener noreferrer">'
        f'<img class="extra-image" src="{escaped_src}" alt="{escaped_value}" '
        f'title="{escaped_value}" loading="lazy" style="max-width: {image_width}px;">'
        f'</a>'
    )


def cell_value(value: str, column: int, link_specs: list[LinkSpec],
               image_columns: set[int], image_width: int) -> str:
    if column in image_columns:
        return image_cell(value, image_width)
    return maybe_link_value(value, column, link_specs)


def as_float(value: str) -> Optional[float]:
    try:
        return float(value)
    except ValueError:
        return None


def color_rule_matches(rule: ColorRule, value: str) -> bool:
    if rule.op in {"lt", "le", "ge", "gt"}:
        lhs = as_float(value)
        rhs = as_float(rule.value)
        if lhs is None or rhs is None:
            return False
        if rule.op == "lt":
            return lhs < rhs
        if rule.op == "le":
            return lhs <= rhs
        if rule.op == "ge":
            return lhs >= rhs
        return lhs > rhs

    if rule.op == "eq":
        return value == rule.value
    if rule.op == "ne":
        return value != rule.value
    if rule.op == "contains":
        return rule.value in value
    if rule.op == "startswith":
        return value.startswith(rule.value)
    if rule.op == "regex":
        return rule.regex.search(value) is not None

    return False


def row_style(extra: list[str], color_column: Optional[int], color_rules: list[ColorRule]) -> str:
    if color_column is None or not color_rules:
        return ""

    value = extra[color_column] if color_column < len(extra) else ""
    for rule in color_rules:
        if color_rule_matches(rule, value):
            return f' style="background-color: {html.escape(rule.color, quote=True)};"'

    return ""


def normalised_rows_for_export(rows, headers):
    result = []
    for row_number, smiles, extra, fields in rows:
        result.append({"row": row_number, "fields": fields})
    return result


def write_export_controls(out, export_filename: Optional[str]):
    if not export_filename:
        return

    out.write('''
<div class="controls">
  <button type="button" id="select-visible">Select visible</button>
  <button type="button" id="clear-selected">Clear selected</button>
  <button type="button" id="download-smiles">Download selected SMILES</button>
  <span id="selected-count">0 selected</span>
</div>
''')


def write_html(rows, display_headers, headers, output, smiles_columns=None, link_specs=None, color_column=None, color_rules=None,
               page_length=50, image_columns=None, image_width=160,
               sanitize=False,
               write_smiles_filename: Optional[str] = None):
    link_specs = link_specs or []
    color_rules = color_rules or []
    image_columns = image_columns or set()
    smiles_columns = smiles_columns or [0]
    smiles_set = set(smiles_columns)

    export_rows_json = json.dumps(normalised_rows_for_export(rows, headers))
    export_filename_json = json.dumps(write_smiles_filename or "selected.smi")

    with open(output, "w") as out:
        out.write("""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>SMILES Browser</title>
<link rel="stylesheet"
 href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
<style>
body {
  font-family: sans-serif;
  margin: 2rem;
}
.controls {
  margin: 1rem 0;
}
.controls button {
  margin-right: 0.5rem;
}
table {
  border-collapse: collapse;
  width: 100%;
}
th, td {
  border: 1px solid #ccc;
  padding: 0.4rem;
  vertical-align: top;
}
th {
  background: #eee;
  position: sticky;
  top: 0;
  z-index: 2;
}
.mol svg {
  max-width: 250px;
  height: auto;
}
.extra-image {
  height: auto;
  display: block;
}
.bad {
  color: red;
  font-weight: bold;
}
.selected-row {
  outline: 3px solid #444;
}
td:first-child,
th:first-child {
  position: sticky;
  left: 0;
  background: #fff;
  z-index: 1;
}
th:first-child {
  background: #eee;
  z-index: 3;
}
</style>
</head>
<body>
<h1>SMILES Browser</h1>
""")

        write_export_controls(out, write_smiles_filename)

        out.write("""<table id="molecules">
<thead>
<tr>
""")

        if write_smiles_filename:
            out.write("<th>Select</th>\n")
        else:
            out.write("<th>#</th>\n")

        for header in display_headers:
            out.write(f"<th>{html.escape(header)}</th>\n")

        out.write("""</tr>
</thead>
<tbody>
""")

        for row_index, (row_number, smiles, extra, fields) in enumerate(rows):
            out.write(f"<tr{row_style(extra, color_column, color_rules)}>\n")
            if write_smiles_filename:
                out.write(
                    f'<td><input type="checkbox" class="row-select" '
                    f'data-row-index="{row_index}" aria-label="Select row {row_number}"></td>\n'
                )
            else:
                out.write(f"<td>{row_number}</td>\n")

            extra_column = 0
            for input_column in range(len(display_headers)):
                value = fields[input_column] if input_column < len(fields) else ""
                if input_column in smiles_set:
                    out.write(f"<td class='mol'>{mol_to_svg(value, sanitize)}</td>\n")
                else:
                    out.write(f"<td>{cell_value(value, extra_column, link_specs, image_columns, image_width)}</td>\n")
                    extra_column += 1

            out.write("</tr>\n")

        out.write("""</tbody>
</table>
<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
<script>
""")
        out.write(f"const moleculeRows = {export_rows_json};\n")
        out.write(f"const smilesDownloadFilename = {export_filename_json};\n")
        out.write(f"const initialPageLength = {page_length};\n")
        out.write(r'''
const selectedRowIndices = new Set();

function setCheckboxState(checkbox, checked) {
    const rowIndex = Number(checkbox.getAttribute('data-row-index'));
    if (Number.isNaN(rowIndex)) {
        return;
    }
    checkbox.checked = checked;
    if (checked) {
        selectedRowIndices.add(rowIndex);
    } else {
        selectedRowIndices.delete(rowIndex);
    }
}

function syncVisibleCheckboxes() {
    document.querySelectorAll('.row-select').forEach(function (checkbox) {
        const rowIndex = Number(checkbox.getAttribute('data-row-index'));
        checkbox.checked = selectedRowIndices.has(rowIndex);
        const row = checkbox.closest('tr');
        if (row) {
            row.classList.toggle('selected-row', checkbox.checked);
        }
    });
}

function selectedIndices() {
    return Array.from(selectedRowIndices).sort(function (a, b) { return a - b; });
}

function updateSelectedCount() {
    $('#selected-count').text(selectedRowIndices.size + ' selected');
}

function smilesLine(row) {
    const smiles = Array.isArray(row.smiles) ? row.smiles : [row.smiles];
    const parts = smiles.concat(row.fields || []);
    return parts.join('\t').replace(/\s+$/g, '');
}

function downloadText(filename, text) {
    const blob = new Blob([text], {type: 'text/plain;charset=utf-8'});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

$(document).ready(function () {
    const table = $('#molecules').DataTable({
        pageLength: initialPageLength
    });

    syncVisibleCheckboxes();
    updateSelectedCount();

    $('#molecules').on('change', '.row-select', function () {
        setCheckboxState(this, this.checked);
        syncVisibleCheckboxes();
        updateSelectedCount();
    });

    // Let users click anywhere in a row, including the molecule image, to toggle selection.
    // Ignore links and form controls so --link-column and --image_column keep working normally.
    $('#molecules tbody').on('click', 'tr', function (event) {
        if ($(event.target).closest('a, button, input, select, textarea').length) {
            return;
        }
        const checkbox = $(this).find('.row-select').get(0);
        if (!checkbox) {
            return;
        }
        setCheckboxState(checkbox, !checkbox.checked);
        syncVisibleCheckboxes();
        updateSelectedCount();
    });

    $('#select-visible').on('click', function () {
        table.rows({search: 'applied', page: 'current'}).nodes().to$()
            .find('.row-select').each(function () {
                setCheckboxState(this, true);
            });
        syncVisibleCheckboxes();
        updateSelectedCount();
    });

    $('#clear-selected').on('click', function () {
        selectedRowIndices.clear();
        syncVisibleCheckboxes();
        updateSelectedCount();
    });

    $('#download-smiles').on('click', function () {
        const indices = selectedIndices();
        if (indices.length === 0) {
            alert('No rows selected');
            return;
        }

        const text = indices.map(i => smilesLine(moleculeRows[i])).join('\n') + '\n';
        downloadText(smilesDownloadFilename, text);
    });

    table.on('draw', function () {
        syncVisibleCheckboxes();
    });
});
</script>
</body>
</html>
''')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input file containing one or more SMILES columns; use - for standard input.")
    parser.add_argument("-o", "--output", default="molecules.html")
    parser.add_argument("-n", type=int, default=None,
                        help="Maximum number of molecules to read")
    parser.add_argument("-s", "--separator", default=None,
                        help="Input separator. Use '\\t' for tab. Defaults to comma for .csv, otherwise whitespace (including stdin).")
    parser.add_argument("--header", action="store_true",
                        help="First record contains column titles.")
    parser.add_argument("--smiles_column", action="append", default=None, metavar="COLUMN[,COLUMN...]",
                        help=("1-based input column(s) containing SMILES. Use a comma-separated list "
                              "or repeat the option. Default is 1."))
    parser.add_argument("--link-column", action="append", default=[], metavar="COLUMN=URL_TEMPLATE",
                        help=("Make a data column a hyperlink. COLUMN is a 1-based field number "
                              "after the structure column, or a header name. The URL template "
                              "may contain {}, {value}, {raw}, or {quoted}. May be repeated."))
    parser.add_argument("--color-column", default=None, metavar="COLUMN",
                        help="Column used for row coloring. Use a 1-based field number or a header name.")
    parser.add_argument("--color-rule", action="append", default=[], metavar="OP:VALUE:COLOR",
                        help=("Row coloring rule. OP is one of lt, le, eq, ne, ge, gt, "
                              "contains, startswith, regex. COLOR is any CSS color. "
                              "Rules are applied in order; first match wins. May be repeated."))
    parser.add_argument("--page-length", type=int, default=50,
                        help="Initial DataTables page length")
    parser.add_argument("--write_smiles", nargs="?", const="selected.smi", default=None,
                        metavar="FILENAME",
                        help=("Add row-selection controls and a browser-side download button. "
                              "The optional FILENAME is the downloaded SMILES filename; "
                              "default is selected.smi."))
    parser.add_argument("--image_column", action="append", default=[], metavar="COLUMN",
                        help=("Render a data column as image thumbnails. COLUMN is a 1-based "
                              "field number after the structure column, or a header name. "
                              "May be repeated."))
    parser.add_argument("--image_width", type=int, default=160,
                        help="Maximum width in pixels for --image_column thumbnails")

    parser.add_argument("--sanitize", action=argparse.BooleanOptionalAction, default=False,
                        help="Sanitize molecules upon read.")
    args = parser.parse_args()

    if args.n is not None and args.n < 1:
        sys.exit("-n must be a positive integer")
    if args.page_length < 1:
        sys.exit("--page-length must be a positive integer")
    if args.image_width < 1:
        sys.exit("--image_width must be a positive integer")

    separator = separator_from_args(args)
    smiles_columns = parse_smiles_columns(args.smiles_column)
    rows, display_headers, headers = read_smiles_file(
        args.input,
        separator=separator,
        has_header=args.header,
        max_molecules=args.n,
        smiles_columns=smiles_columns,
    )

    link_specs = parse_link_columns(args.link_column, headers)
    image_columns = parse_image_columns(args.image_column, headers)

    color_column = None
    color_rules = []
    if args.color_column is not None:
        color_column = column_from_string(args.color_column, headers, "--color-column")
        color_rules = parse_color_rules(args.color_rule)
        if not color_rules:
            sys.exit("--color-column requires at least one --color-rule")
    elif args.color_rule:
        sys.exit("--color-rule requires --color-column")

    write_html(
        rows,
        display_headers,
        headers,
        args.output,
        smiles_columns=smiles_columns,
        link_specs=link_specs,
        color_column=color_column,
        color_rules=color_rules,
        page_length=args.page_length,
        image_columns=image_columns,
        image_width=args.image_width,
        sanitize=args.sanitize,
        write_smiles_filename=args.write_smiles,
    )


if __name__ == "__main__":
    main()
