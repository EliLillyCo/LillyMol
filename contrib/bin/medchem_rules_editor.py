#!/usr/bin/env python3
"""Generate a static HTML editor for Lilly Medchem Rules demerit values.

The generated page is self-contained. It lets a user review each active query,
set a replacement demerit value, ignore a rule by assigning zero, turn
a rule into a hard reject by assigning a configured high value, and download a
complete query-name/value override file.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import html
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError as err:  # pragma: no cover - depends on local environment.
    Chem = None
    rdDepictor = None
    rdMolDraw2D = None
    RDKIT_IMPORT_ERROR = err
else:
    RDKIT_IMPORT_ERROR = None


DEFAULT_QUERY_DIR = Path("data/queries/LillyMedchemRules")
DEFAULT_MANIFEST = "demerits"
DEFAULT_OUTPUT = Path("medchem_rules_editor.html")
DEFAULT_EXAMPLES = 3
DEFAULT_MAX_DEMERIT = 100
DEFAULT_REJECT_VALUE = 100
DEFAULT_IMAGE_WIDTH = 220
DEFAULT_IMAGE_HEIGHT = 150


@dataclass
class Example:
    smiles: str
    name: str
    svg: str | None
    error: str | None


@dataclass
class Rule:
    index: int
    filename: str
    query_name: str
    original_demerit: int
    examples: list[Example]


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a static HTML editor for Lilly Medchem Rules demerits."
    )
    parser.add_argument(
        "--query-dir",
        type=Path,
        default=DEFAULT_QUERY_DIR,
        help=f"Directory containing query, .smi, and manifest files. Default: {DEFAULT_QUERY_DIR}",
    )
    parser.add_argument(
        "--manifest",
        default=DEFAULT_MANIFEST,
        help=f"Manifest file name or path. Default: {DEFAULT_MANIFEST}",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output HTML file. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--examples",
        type=int,
        default=DEFAULT_EXAMPLES,
        help=f"Maximum example molecules to show per query. Default: {DEFAULT_EXAMPLES}",
    )
    parser.add_argument(
        "--max-demerit",
        type=int,
        default=DEFAULT_MAX_DEMERIT,
        help=f"Maximum allowed editable demerit value. Default: {DEFAULT_MAX_DEMERIT}",
    )
    parser.add_argument(
        "--reject-value",
        type=int,
        default=DEFAULT_REJECT_VALUE,
        help=f"Value assigned by the Hard reject button. Default: {DEFAULT_REJECT_VALUE}",
    )
    parser.add_argument(
        "--title",
        default="Lilly Medchem Rules Demerit Editor",
        help="Title shown in the generated HTML page.",
    )
    parser.add_argument(
        "--image-width",
        type=int,
        default=DEFAULT_IMAGE_WIDTH,
        help=f"SVG depiction width in pixels. Default: {DEFAULT_IMAGE_WIDTH}",
    )
    parser.add_argument(
        "--image-height",
        type=int,
        default=DEFAULT_IMAGE_HEIGHT,
        help=f"SVG depiction height in pixels. Default: {DEFAULT_IMAGE_HEIGHT}",
    )
    parser.add_argument(
        "--require-rdkit",
        action="store_true",
        help="Fail if RDKit is not available instead of falling back to SMILES text.",
    )
    parser.add_argument(
        "--strict-filenames",
        action="store_true",
        help=(
            "Fail if a query name does not match the query filename stem. "
            "Hyphens and underscores are treated as different."
        ),
    )
    return parser.parse_args(argv)


def manifest_path(query_dir: Path, manifest: str) -> Path:
    path = Path(manifest)
    if path.is_absolute() or len(path.parts) > 1:
        return path
    return query_dir / path


def read_manifest(path: Path) -> list[str]:
    result: list[str] = []
    with path.open("r", encoding="utf-8") as reader:
        for line_number, line in enumerate(reader, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "#" in line:
                line = line.split("#", 1)[0].strip()
            if not line:
                continue
            name = Path(line.split()[0]).name
            if not name.endswith((".qry", ".textproto")):
                raise ValueError(f"{path}:{line_number}: expected .qry or .textproto entry, got {name!r}")
            result.append(name)
    return result


def parse_query_file(path: Path) -> tuple[str, int]:
    text = path.read_text(encoding="utf-8", errors="replace")

    name_match = re.search(r'^\s*name:\s*"([^"]+)"', text, re.MULTILINE)
    value_match = re.search(r'^\s*numeric_value:\s*(\d+)', text, re.MULTILINE)
    if name_match or value_match:
        if not name_match:
            raise ValueError(f"{path}: textproto query has no name field")
        if not value_match:
            raise ValueError(f"{path}: textproto query has no numeric_value field")
        return name_match.group(1), int(value_match.group(1))

    name_match = re.search(r'\(A\s+C\s+comment\s+"([^"]+)"\)', text, re.IGNORECASE)
    value_match = re.search(r'\(A\s+D\s+numeric_value\s+(\d+)\)', text, re.IGNORECASE)
    if not name_match:
        raise ValueError(f"{path}: old-style query has no comment field")
    if not value_match:
        raise ValueError(f"{path}: old-style query has no numeric_value field")
    return name_match.group(1), int(value_match.group(1))


def mol_to_svg(smiles: str, width: int, height: int, require_rdkit: bool) -> tuple[str | None, str | None]:
    if RDKIT_IMPORT_ERROR is not None:
        if require_rdkit:
            raise RuntimeError("RDKit is required to generate molecule depictions") from RDKIT_IMPORT_ERROR
        return None, "RDKit not available"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES"

    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = svg.replace("<?xml version='1.0' encoding='iso-8859-1'?>", "")
    return svg, None


def read_examples(path: Path, limit: int, width: int, height: int, require_rdkit: bool) -> list[Example]:
    if limit <= 0:
        return []

    examples: list[Example] = []
    with path.open("r", encoding="utf-8", errors="replace") as reader:
        for line in reader:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split(maxsplit=1)
            smiles = tokens[0]
            name = tokens[1] if len(tokens) > 1 else ""
            svg, error = mol_to_svg(smiles, width, height, require_rdkit)
            examples.append(Example(smiles=smiles, name=name, svg=svg, error=error))
            if len(examples) >= limit:
                break
    return examples


def build_rules(args: argparse.Namespace) -> tuple[list[Rule], list[str]]:
    qdir = args.query_dir
    manifest = manifest_path(qdir, args.manifest)
    entries = read_manifest(manifest)
    rules: list[Rule] = []
    warnings: list[str] = []
    seen_names: dict[str, str] = {}

    for index, filename in enumerate(entries, start=1):
        query_path = qdir / filename
        if not query_path.is_file():
            raise FileNotFoundError(f"Query file listed in {manifest} not found: {query_path}")

        query_name, demerit = parse_query_file(query_path)
        stem = query_path.stem
        if query_name != stem:
            message = f"{filename}: query name {query_name!r} does not match filename stem {stem!r}"
            if args.strict_filenames:
                raise ValueError(message)
            warnings.append(message)

        if query_name in seen_names:
            raise ValueError(
                f"Duplicate query name {query_name!r}: {seen_names[query_name]} and {filename}"
            )
        seen_names[query_name] = filename

        if demerit < 0 or demerit > args.max_demerit:
            warnings.append(
                f"{filename}: original demerit {demerit} is outside 0..{args.max_demerit}"
            )

        smi_path = query_path.with_suffix(".smi")
        if not smi_path.is_file():
            raise FileNotFoundError(f"Expected example file not found: {smi_path}")

        rules.append(
            Rule(
                index=index,
                filename=filename,
                query_name=query_name,
                original_demerit=demerit,
                examples=read_examples(
                    smi_path,
                    args.examples,
                    args.image_width,
                    args.image_height,
                    args.require_rdkit,
                ),
            )
        )

    return rules, warnings


def example_to_html(example: Example) -> str:
    parts: list[str] = ['<div class="example">']
    if example.svg:
        parts.append(f'<div class="mol-svg">{example.svg}</div>')
    else:
        parts.append('<div class="mol-fallback">')
        parts.append(f'<code>{html.escape(example.smiles)}</code>')
        if example.error:
            parts.append(f'<span>{html.escape(example.error)}</span>')
        parts.append('</div>')

    parts.append('<div class="example-meta">')
    parts.append(f'<code>{html.escape(example.smiles)}</code>')
    if example.name:
        parts.append(f'<span>{html.escape(example.name)}</span>')
    parts.append('</div></div>')
    return "".join(parts)


def rules_for_json(rules: list[Rule]) -> list[dict[str, Any]]:
    return [
        {
            "index": rule.index,
            "filename": rule.filename,
            "name": rule.query_name,
            "original": rule.original_demerit,
        }
        for rule in rules
    ]


def build_html(args: argparse.Namespace, rules: list[Rule], warnings: list[str]) -> str:
    now = _datetime.datetime.now().isoformat(timespec="seconds")
    qdir = str(args.query_dir)
    manifest = str(manifest_path(args.query_dir, args.manifest))
    rules_json = json.dumps(rules_for_json(rules), ensure_ascii=False)
    warnings_json = json.dumps(warnings, ensure_ascii=False)
    title = html.escape(args.title)

    parts: list[str] = []
    parts.append("""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
""")
    parts.append(f"<title>{title}</title>\n")
    parts.append("""<style>
:root {
  --bg: #f6f7f9;
  --panel: #ffffff;
  --ink: #20242a;
  --muted: #69707d;
  --border: #d7dce3;
  --changed: #fff5cc;
  --ignore: #edf7ed;
  --hard-reject: #fde8e8;
  --accent: #246b73;
  --accent-dark: #16484e;
  --danger: #a83232;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: Arial, Helvetica, sans-serif;
  background: var(--bg);
  color: var(--ink);
  font-size: 14px;
}
header {
  padding: 18px 24px 12px 24px;
  background: var(--panel);
  border-bottom: 1px solid var(--border);
}
h1 {
  margin: 0 0 6px 0;
  font-size: 22px;
  line-height: 1.2;
}
.summary {
  margin: 0;
  color: var(--muted);
  line-height: 1.4;
}
.toolbar {
  position: sticky;
  top: 0;
  z-index: 20;
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 10px;
  padding: 10px 24px;
  background: #edf1f4;
  border-bottom: 1px solid var(--border);
}
.toolbar input[type="search"] {
  width: min(420px, 100%);
  padding: 7px 9px;
  border: 1px solid var(--border);
  border-radius: 4px;
  font-size: 14px;
}
.toolbar label {
  display: inline-flex;
  align-items: center;
  gap: 5px;
  white-space: nowrap;
  color: #343a42;
}
button {
  border: 1px solid var(--border);
  border-radius: 4px;
  background: #ffffff;
  color: #20242a;
  padding: 7px 10px;
  font: inherit;
  cursor: pointer;
}
button:hover { border-color: #9ca5b1; }
button.primary {
  background: var(--accent);
  border-color: var(--accent);
  color: #fff;
}
button.primary:hover { background: var(--accent-dark); }
button.ignore {
  border-color: #9fc5a1;
  color: #21612a;
}
button.hard-reject {
  border-color: #e1aaaa;
  color: var(--danger);
}
.status {
  margin-left: auto;
  color: var(--muted);
  white-space: nowrap;
}
main { padding: 16px 24px 32px 24px; }
.warning-box {
  margin-bottom: 14px;
  padding: 10px 12px;
  border: 1px solid #dbb75b;
  background: #fff8df;
}
.warning-box[hidden] { display: none; }
.table-wrap {
  overflow-x: auto;
  border: 1px solid var(--border);
  background: var(--panel);
}
table {
  width: 100%;
  border-collapse: collapse;
  min-width: 1050px;
}
th, td {
  border-bottom: 1px solid var(--border);
  padding: 8px 10px;
  vertical-align: top;
  text-align: left;
}
th {
  background: #dfe6eb;
  font-size: 12px;
  text-transform: uppercase;
  color: #3d4651;
  letter-spacing: 0;
}
tbody tr.changed { background: var(--changed); }
tbody tr.ignored { background: var(--ignore); }
tbody tr.hard-rejected { background: var(--hard-reject); }
.index { width: 52px; color: var(--muted); text-align: right; }
.name-cell strong { display: block; font-size: 14px; }
.name-cell code {
  display: block;
  margin-top: 3px;
  color: var(--muted);
  font-size: 12px;
  overflow-wrap: anywhere;
}
.original-cell { width: 95px; text-align: right; font-variant-numeric: tabular-nums; }
.action-cell { width: 330px; white-space: nowrap; }
.demerit-controls {
  display: inline-flex;
  align-items: center;
  gap: 6px;
}
.demerit-controls input {
  width: 72px;
  padding: 6px 7px;
  border: 1px solid var(--border);
  border-radius: 4px;
  font: inherit;
}
.demerit-controls input.invalid {
  border-color: var(--danger);
  outline: 1px solid var(--danger);
}
.examples-cell { min-width: 520px; }
.examples {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
}
.example {
  width: 220px;
  border: 1px solid #e2e5e9;
  background: #fff;
  padding: 6px;
}
.mol-svg {
  width: 100%;
  min-height: 150px;
  display: flex;
  align-items: center;
  justify-content: center;
  overflow: hidden;
}
.mol-svg svg { max-width: 100%; height: auto; }
.mol-fallback {
  min-height: 150px;
  display: flex;
  flex-direction: column;
  justify-content: center;
  gap: 8px;
  color: var(--muted);
  overflow-wrap: anywhere;
}
.example-meta {
  border-top: 1px solid #edf0f2;
  padding-top: 5px;
  color: var(--muted);
  font-size: 11px;
  line-height: 1.35;
}
.example-meta code {
  display: block;
  color: #27313a;
  overflow-wrap: anywhere;
}
.example-meta span {
  display: block;
  overflow-wrap: anywhere;
}
.no-examples { color: var(--muted); }
.hidden-row { display: none; }
@media (max-width: 720px) {
  header, .toolbar, main { padding-left: 12px; padding-right: 12px; }
  .status { width: 100%; margin-left: 0; }
}
</style>
</head>
<body>
""")
    parts.append("<header>\n")
    parts.append(f"<h1>{title}</h1>\n")
    parts.append(
        '<p class="summary">'
        f'{len(rules)} rules from <code>{html.escape(manifest)}</code>. '
        f'Examples per rule: {args.examples}. Max demerit: {args.max_demerit}. '
        f'Hard reject value: {args.reject_value}. Generated: {html.escape(now)}.'
        '</p>\n'
    )
    parts.append("</header>\n")
    parts.append("""<div class="toolbar">
<input id="search" type="search" placeholder="Filter by query name, filename, SMILES, or molecule id" oninput="updateView()">
<label><input id="changedOnly" type="checkbox" onchange="updateView()"> Changed</label>
<label><input id="ignoredOnly" type="checkbox" onchange="updateView()"> Ignored</label>
<label><input id="hardRejectedOnly" type="checkbox" onchange="updateView()"> Hard reject</label>
<button type="button" class="primary" onclick="downloadOverrides()">Download</button>
<button type="button" onclick="copyOverrides()">Copy</button>
<button type="button" onclick="resetVisible()">Reset visible</button>
<button type="button" onclick="resetAll()">Reset all</button>
<span id="status" class="status"></span>
</div>
<main>
<div id="warningBox" class="warning-box" hidden></div>
<div class="table-wrap">
<table>
<thead>
<tr>
<th class="index">#</th>
<th>Query</th>
<th class="original-cell">Original</th>
<th class="action-cell">Demerit</th>
<th class="examples-cell">Examples</th>
</tr>
</thead>
<tbody>
""")

    for rule in rules:
        searchable = " ".join(
            [rule.query_name, rule.filename]
            + [example.smiles for example in rule.examples]
            + [example.name for example in rule.examples]
        ).lower()
        parts.append(
            f'<tr data-index="{rule.index}" data-name="{html.escape(rule.query_name, quote=True)}" '
            f'data-original="{rule.original_demerit}" data-search="{html.escape(searchable, quote=True)}">\n'
        )
        parts.append(f'<td class="index">{rule.index}</td>\n')
        parts.append('<td class="name-cell">')
        parts.append(f'<strong>{html.escape(rule.query_name)}</strong>')
        parts.append(f'<code>{html.escape(rule.filename)}</code>')
        parts.append('</td>\n')
        parts.append(f'<td class="original-cell">{rule.original_demerit}</td>\n')
        parts.append(
            '<td class="action-cell"><div class="demerit-controls">'
            f'<button type="button" class="ignore" onclick="ignoreRule({rule.index})">Ignore</button> '
            f'<input type="number" min="0" max="{args.max_demerit}" step="1" '
            f'value="{rule.original_demerit}" data-index="{rule.index}" oninput="updateView()"> '
            f'<button type="button" class="hard-reject" onclick="hardRejectRule({rule.index})">Hard reject</button> '
            f'<button type="button" onclick="resetRule({rule.index})">Reset</button>'
            '</div></td>\n'
        )
        parts.append('<td class="examples-cell">')
        if rule.examples:
            parts.append('<div class="examples">')
            for example in rule.examples:
                parts.append(example_to_html(example))
            parts.append('</div>')
        else:
            parts.append('<span class="no-examples">No examples</span>')
        parts.append('</td>\n</tr>\n')

    parts.append("""</tbody>
</table>
</div>
</main>
<script>
""")
    parts.append(f"const RULES = {rules_json};\n")
    parts.append(f"const WARNINGS = {warnings_json};\n")
    parts.append(f"const MAX_DEMERIT = {args.max_demerit};\n")
    parts.append(f"const REJECT_VALUE = {args.reject_value};\n")
    parts.append(f"const SOURCE_MANIFEST = {json.dumps(manifest)};\n")
    parts.append(f"const SOURCE_QUERY_DIR = {json.dumps(qdir)};\n")
    parts.append(r'''
function inputFor(index) {
  return document.querySelector(`input[data-index="${index}"]`);
}

function rowFor(index) {
  return document.querySelector(`tr[data-index="${index}"]`);
}

function valueFor(index) {
  const input = inputFor(index);
  return Number.parseInt(input.value, 10);
}

function isValidInput(input) {
  if (!/^\d+$/.test(input.value)) {
    return false;
  }
  const value = Number.parseInt(input.value, 10);
  return value >= 0 && value <= MAX_DEMERIT;
}

function isChanged(rule) {
  return valueFor(rule.index) !== rule.original;
}

function isIgnored(rule) {
  return valueFor(rule.index) === 0;
}

function isHardRejected(rule) {
  return valueFor(rule.index) === REJECT_VALUE;
}

function updateView() {
  const query = document.getElementById("search").value.trim().toLowerCase();
  const changedOnly = document.getElementById("changedOnly").checked;
  const ignoredOnly = document.getElementById("ignoredOnly").checked;
  const hardRejectedOnly = document.getElementById("hardRejectedOnly").checked;
  let changed = 0;
  let ignored = 0;
  let hardRejected = 0;
  let invalid = 0;
  let visible = 0;

  for (const rule of RULES) {
    const row = rowFor(rule.index);
    const input = inputFor(rule.index);
    const valid = isValidInput(input);
    input.classList.toggle("invalid", !valid);
    if (!valid) {
      invalid += 1;
    }

    const changedRow = valid && isChanged(rule);
    const ignoredRow = valid && isIgnored(rule);
    const hardRejectedRow = valid && isHardRejected(rule);
    if (changedRow) {
      changed += 1;
    }
    if (ignoredRow) {
      ignored += 1;
    }
    if (hardRejectedRow) {
      hardRejected += 1;
    }

    row.classList.toggle("changed", changedRow);
    row.classList.toggle("ignored", ignoredRow);
    row.classList.toggle("hard-rejected", hardRejectedRow);

    const searchText = row.getAttribute("data-search") || "";
    const matchesSearch = query === "" || searchText.includes(query);
    const matchesChanged = !changedOnly || changedRow;
    const matchesIgnored = !ignoredOnly || ignoredRow;
    const matchesHardRejected = !hardRejectedOnly || hardRejectedRow;
    const show = matchesSearch && matchesChanged && matchesIgnored && matchesHardRejected;
    row.classList.toggle("hidden-row", !show);
    if (show) {
      visible += 1;
    }
  }

  const status = document.getElementById("status");
  status.textContent = `${visible} visible; ${changed} changed; ${ignored} ignored; ${hardRejected} hard reject; ${invalid} invalid`;

  for (const button of document.querySelectorAll("button.primary")) {
    button.disabled = invalid > 0;
  }
}

function ignoreRule(index) {
  inputFor(index).value = 0;
  updateView();
}

function hardRejectRule(index) {
  inputFor(index).value = REJECT_VALUE;
  updateView();
}

function resetRule(index) {
  const rule = RULES[index - 1];
  inputFor(index).value = rule.original;
  updateView();
}

function resetVisible() {
  for (const rule of RULES) {
    const row = rowFor(rule.index);
    if (!row.classList.contains("hidden-row")) {
      inputFor(rule.index).value = rule.original;
    }
  }
  updateView();
}

function resetAll() {
  for (const rule of RULES) {
    inputFor(rule.index).value = rule.original;
  }
  updateView();
}

function buildOverrides() {
  const lines = [
    `# Generated by medchem_rules_editor.py`,
    `# source_query_dir: ${SOURCE_QUERY_DIR}`,
    `# source_manifest: ${SOURCE_MANIFEST}`,
    `# max_demerit: ${MAX_DEMERIT}`,
    `# reject_value: ${REJECT_VALUE}`,
  ];
  for (const rule of RULES) {
    const input = inputFor(rule.index);
    if (!isValidInput(input)) {
      throw new Error(`Invalid value for ${rule.name}`);
    }
    lines.push(`${rule.name} ${Number.parseInt(input.value, 10)}`);
  }
  return lines.join("\n") + "\n";
}

function downloadOverrides() {
  const text = buildOverrides();
  const blob = new Blob([text], {type: "text/plain;charset=utf-8"});
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "medchem_rules_demerits.txt";
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}

async function copyOverrides() {
  const text = buildOverrides();
  await navigator.clipboard.writeText(text);
}

function showWarnings() {
  if (WARNINGS.length === 0) {
    return;
  }
  const box = document.getElementById("warningBox");
  box.hidden = false;
  box.innerHTML = `<strong>${WARNINGS.length} generation warning(s)</strong><br>` +
    WARNINGS.map((warning) => `<div>${warning.replace(/[&<>"]/g, (ch) => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;"}[ch]))}</div>`).join("");
}

showWarnings();
updateView();
''')
    parts.append("</script>\n</body>\n</html>\n")
    return "".join(parts)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if args.examples < 0:
        raise ValueError("--examples must be non-negative")
    if args.max_demerit < 1:
        raise ValueError("--max-demerit must be positive")
    if args.reject_value < 0 or args.reject_value > args.max_demerit:
        raise ValueError("--reject-value must be between 0 and --max-demerit")

    rules, warnings = build_rules(args)
    page = build_html(args, rules, warnings)
    args.output.write_text(page, encoding="utf-8")
    print(f"Wrote {args.output} with {len(rules)} rules", file=sys.stderr)
    if warnings:
        print(f"Generated with {len(warnings)} warning(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
