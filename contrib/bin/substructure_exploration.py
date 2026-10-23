#!/usr/bin/env python3
"""Explore activity distributions associated with substructure motifs."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from google.protobuf import text_format

from substructure_exploration_pb2 import SubstructureExploration

from lillymol import Molecule
from lillymol_io import ReaderContext
from lillymol_query import SubstructureQuery


@dataclass
class MoleculeRecord:
    smiles: str
    identifier: str
    activity: float
    mol: Molecule


@dataclass
class Example:
    smiles: str
    labelled_smiles: str
    identifier: str
    activity: float
    nhits: int


@dataclass
class QuerySpec:
    name: str
    source: str
    source_type: str
    query: SubstructureQuery
    ncon: list[int] = field(default_factory=list)
    charge_assigner: bool = False


@dataclass
class QueryAccumulator:
    spec: QuerySpec
    matched: list[float] = field(default_factory=list)
    unmatched: list[float] = field(default_factory=list)
    count_values: dict[str, list[float]] = field(default_factory=dict)
    examples: list[Example] = field(default_factory=list)
    molecules_with_multiple_matches: int = 0
    max_matches: int = 0


def die(message: str) -> None:
    print(message, file=sys.stderr)
    sys.exit(1)


def separator_for_file(path: str) -> str | None:
    if path.lower().endswith('.csv'):
        return ','
    return None


def read_config(path: str) -> SubstructureExploration:
    message = SubstructureExploration()
    try:
        text = Path(path).read_text()
    except OSError as exc:
        die(f"Cannot read config '{path}': {exc}")

    try:
        text_format.Parse(text, message)
    except text_format.ParseError as exc:
        die(f"Cannot parse config '{path}': {exc}")

    return message


def build_queries(config: SubstructureExploration) -> list[QuerySpec]:
    result = []
    for ndx, proto in enumerate(config.query, 1):
        qry = SubstructureQuery()
        source_type = proto.WhichOneof('query_spec')
        if source_type == 'smarts':
            source = proto.smarts
            if not qry.build_from_smarts(source):
                die(f"Query {ndx}: invalid SMARTS '{source}'")
        elif source_type == 'fname':
            source = proto.fname
            if source.endswith('.qry'):
                ok = qry.read_msi(source)
            else:
                ok = qry.read_proto(source) or qry.read_msi(source)
            if not ok:
                die(f"Query {ndx}: cannot read query file '{source}'")
        else:
            die(f"Query {ndx}: no smarts or fname specified")

        if not proto.HasField('name') or not proto.name:
            die(f"Query {ndx}: every query must have a name")

        result.append(QuerySpec(
            name=proto.name,
            source=source,
            source_type=source_type,
            query=qry,
            ncon=list(proto.ncon),
            charge_assigner=proto.charge_assigner if proto.HasField('charge_assigner') else False,
        ))

    names = [q.name for q in result]
    duplicate_names = sorted(name for name in set(names) if names.count(name) > 1)
    if duplicate_names:
        die(f"Duplicate query names: {', '.join(duplicate_names)}")

    return result


def read_activity(path: str) -> dict[str, float]:
    result: dict[str, float] = {}
    sep = separator_for_file(path)
    try:
        with open(path, newline='') as reader:
            if sep == ',':
                records = csv.reader(reader)
                header = next(records, None)
                if header is None:
                    die(f"{path}: empty activity file")
                for line_number, row in enumerate(records, 2):
                    if not row:
                        continue
                    if len(row) < 2:
                        die(f"{path}:{line_number}: activity records need at least two columns")
                    add_activity_record(result, row[0], row[1], path, line_number)
            else:
                header = reader.readline()
                if not header:
                    die(f"{path}: empty activity file")
                for line_number, line in enumerate(reader, 2):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    fields = line.split()
                    if len(fields) < 2:
                        die(f"{path}:{line_number}: activity records need at least two columns")
                    add_activity_record(result, fields[0], fields[1], path, line_number)
    except OSError as exc:
        die(f"Cannot read activity file '{path}': {exc}")

    if not result:
        die(f"{path}: no activity values read")

    return result


def add_activity_record(activity: dict[str, float], identifier: str, value: str,
                        fname: str, line_number: int) -> None:
    identifier = identifier.strip()
    if not identifier:
        die(f"{fname}:{line_number}: empty identifier")
    if identifier in activity:
        die(f"{fname}:{line_number}: duplicate activity for '{identifier}'")
    try:
        activity[identifier] = float(value)
    except ValueError:
        die(f"{fname}:{line_number}: invalid activity '{value}'")


def read_molecules(path: str, activity: dict[str, float], okmissing_activity: bool) -> tuple[list[MoleculeRecord], int]:
    result: list[MoleculeRecord] = []
    missing_activity = 0
    seen: set[str] = set()

    with ReaderContext(path) as reader:
        for mol in reader:
            identifier = str(mol.name()).split()[0]
            if not identifier:
                die(f"{path}: molecule with no identifier")
            if identifier in seen:
                die(f"{path}: duplicate molecule identifier '{identifier}'")
            seen.add(identifier)

            if identifier not in activity:
                if not okmissing_activity:
                    die(f"{path}: missing activity for '{identifier}'")
                missing_activity += 1
                continue

            result.append(MoleculeRecord(
                smiles=str(mol.smiles()),
                identifier=identifier,
                activity=activity[identifier],
                mol=mol,
            ))

    if not result:
        die(f"{path}: no molecules with activity")

    unused_activity = sorted(set(activity) - seen)
    if unused_activity:
        print(f"Warning: {len(unused_activity)} activity records not found in smiles file", file=sys.stderr)

    return result, missing_activity


def percentile(sorted_values: list[float], pct: float) -> float | None:
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = (len(sorted_values) - 1) * pct / 100.0
    lower = math.floor(pos)
    upper = math.ceil(pos)
    if lower == upper:
        return sorted_values[int(pos)]
    fraction = pos - lower
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction


def summary(values: Iterable[float]) -> dict[str, float | int | None]:
    vals = sorted(values)
    if not vals:
        return {
            'n': 0,
            'min': None,
            'q1': None,
            'median': None,
            'mean': None,
            'q3': None,
            'max': None,
        }

    return {
        'n': len(vals),
        'min': vals[0],
        'q1': percentile(vals, 25),
        'median': percentile(vals, 50),
        'mean': statistics.fmean(vals),
        'q3': percentile(vals, 75),
        'max': vals[-1],
    }


def common_bins(values: list[float], nbins: int) -> dict[str, object]:
    minval = min(values)
    maxval = max(values)
    if minval == maxval:
        edges = [minval, maxval]
        return {'min': minval, 'max': maxval, 'edges': edges}

    nbins = max(1, nbins)
    width = (maxval - minval) / nbins
    edges = [minval + i * width for i in range(nbins)] + [maxval]
    return {'min': minval, 'max': maxval, 'edges': edges}


def histogram(values: Iterable[float], edges: list[float]) -> list[int]:
    counts = [0] * (len(edges) - 1)
    if not counts:
        return counts

    minval = edges[0]
    maxval = edges[-1]
    if minval == maxval:
        counts[0] = sum(1 for _ in values)
        return counts

    width = (maxval - minval) / len(counts)
    for value in values:
        bucket = int((value - minval) / width)
        if bucket < 0:
            bucket = 0
        elif bucket >= len(counts):
            bucket = len(counts) - 1
        counts[bucket] += 1
    return counts


def bucket_for_count(nhits: int, buckets: list[str]) -> str:
    for bucket in buckets:
        if bucket.endswith('+'):
            try:
                threshold = int(bucket[:-1])
            except ValueError:
                continue
            if nhits >= threshold:
                return bucket
        else:
            try:
                if nhits == int(bucket):
                    return bucket
            except ValueError:
                continue
    return buckets[-1]


def parse_count_buckets(text: str) -> list[str]:
    buckets = [token.strip() for token in text.split(',') if token.strip()]
    if not buckets:
        die('--count-buckets must specify at least one bucket')
    if buckets[0] != '0':
        die('--count-buckets must include 0 as the first bucket')
    return buckets


def labelled_example(record: MoleculeRecord, matches, nhits: int) -> Example:
    mol = Molecule(record.mol)
    for match in matches or []:
        mol.set_isotopes(match, 1)
    return Example(
        smiles=record.smiles,
        labelled_smiles=str(mol.aromatic_smiles()),
        identifier=record.identifier,
        activity=record.activity,
        nhits=nhits,
    )


def analyse(args: argparse.Namespace) -> dict[str, object]:
    config = read_config(args.config)
    queries = build_queries(config)
    if not queries:
        die(f"{args.config}: no queries")

    unsupported = [q.name for q in queries if q.charge_assigner or q.ncon]
    if unsupported:
        print(
            'Warning: charge_assigner and ncon query metadata are recorded but not applied by this script: ' +
            ', '.join(unsupported),
            file=sys.stderr,
        )

    activity = read_activity(args.activity)
    molecules, missing_activity = read_molecules(args.smiles, activity, args.okmissing_activity)
    bins = common_bins([m.activity for m in molecules], args.bins)
    edges = bins['edges']
    count_buckets = parse_count_buckets(args.count_buckets)

    accumulators = [QueryAccumulator(q) for q in queries]
    for acc in accumulators:
        acc.count_values = {bucket: [] for bucket in count_buckets}

    for record in molecules:
        for acc in accumulators:
            matches = acc.spec.query.substructure_search_matches(record.mol)
            nhits = len(matches) if matches is not None else 0
            acc.max_matches = max(acc.max_matches, nhits)
            if nhits > 1:
                acc.molecules_with_multiple_matches += 1

            bucket = bucket_for_count(nhits, count_buckets)
            acc.count_values[bucket].append(record.activity)

            if nhits > 0:
                acc.matched.append(record.activity)
                if len(acc.examples) < args.max_examples:
                    acc.examples.append(labelled_example(record, matches, nhits))
            else:
                acc.unmatched.append(record.activity)

    results = []
    for ndx, acc in enumerate(accumulators, 1):
        matched_summary = summary(acc.matched)
        unmatched_summary = summary(acc.unmatched)
        delta_mean = None
        if matched_summary['mean'] is not None and unmatched_summary['mean'] is not None:
            delta_mean = matched_summary['mean'] - unmatched_summary['mean']
        delta_median = None
        if matched_summary['median'] is not None and unmatched_summary['median'] is not None:
            delta_median = matched_summary['median'] - unmatched_summary['median']

        results.append({
            'index': ndx,
            'name': acc.spec.name,
            'source_type': acc.spec.source_type,
            'source': acc.spec.source,
            'ncon': acc.spec.ncon,
            'charge_assigner': acc.spec.charge_assigner,
            'matched': {
                'summary': matched_summary,
                'histogram': histogram(acc.matched, edges),
            },
            'unmatched': {
                'summary': unmatched_summary,
                'histogram': histogram(acc.unmatched, edges),
            },
            'delta_mean': delta_mean,
            'delta_median': delta_median,
            'multiplicity': {
                'molecules_with_multiple_matches': acc.molecules_with_multiple_matches,
                'max_matches': acc.max_matches,
                'profiles': {
                    bucket: {
                        'summary': summary(values),
                        'histogram': histogram(values, edges),
                    }
                    for bucket, values in acc.count_values.items()
                },
            },
            'examples': [example.__dict__ for example in acc.examples],
        })

    return {
        'metadata': {
            'smiles': args.smiles,
            'activity': args.activity,
            'config': args.config,
            'molecules_read_with_activity': len(molecules),
            'molecules_missing_activity': missing_activity,
            'queries': len(queries),
            'histogram_bins': args.bins,
            'multiplicity_mode': args.multiplicity,
            'count_buckets': count_buckets,
        },
        'global': {
            'summary': summary(m.activity for m in molecules),
            'histogram': histogram((m.activity for m in molecules), edges),
            'bins': bins,
        },
        'queries': results,
    }


def fmt(value, digits: int = 3) -> str:
    if value is None:
        return ''
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}".rstrip('0').rstrip('.')


def simple_histogram_svg(counts: list[int], width: int = 220, height: int = 44) -> str:
    max_count = max(counts + [1])
    nbins = max(len(counts), 1)
    slot = width / nbins
    bar_width = max(1.0, slot * 0.72)
    parts = [f'<svg viewBox="0 0 {width} {height}" class="hist global-hist" aria-hidden="true">']
    parts.append(f'<line x1="0" y1="{height - 1}" x2="{width}" y2="{height - 1}"/>')
    for i, count in enumerate(counts):
        h = (count / max_count) * (height - 6)
        x = i * slot + slot * 0.14
        y = height - h - 1
        parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{h:.2f}"/>')
    parts.append('</svg>')
    return ''.join(parts)


def histogram_svg(matched: list[int], unmatched: list[int], width: int = 220, height: int = 54) -> str:
    max_count = max(matched + unmatched + [1])
    nbins = max(len(matched), 1)
    slot = width / nbins
    bar_width = max(1.0, slot * 0.38)
    parts = [f'<svg viewBox="0 0 {width} {height}" class="hist" aria-hidden="true">']
    parts.append(f'<line x1="0" y1="{height - 1}" x2="{width}" y2="{height - 1}"/>')
    for i, count in enumerate(unmatched):
        h = (count / max_count) * (height - 6)
        x = i * slot + slot * 0.12
        y = height - h - 1
        parts.append(f'<rect class="unmatched" x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{h:.2f}"/>')
    for i, count in enumerate(matched):
        h = (count / max_count) * (height - 6)
        x = i * slot + slot * 0.50
        y = height - h - 1
        parts.append(f'<rect class="matched" x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{h:.2f}"/>')
    parts.append('</svg>')
    return ''.join(parts)


def multiplicity_profiles_html(row: dict[str, object]) -> str:
    profiles = row['multiplicity']['profiles']
    parts = ['<div class="count-profiles">']
    for bucket, profile in profiles.items():
        s = profile['summary']
        if s['n'] == 0:
            continue
        parts.append(
            f'<span><b>{html.escape(bucket)}</b> n {s["n"]} med {fmt(s["median"])}</span>'
        )
    parts.append('</div>')
    return ''.join(parts)


def examples_html(examples: list[dict[str, object]]) -> str:
    if not examples:
        return '<span class="muted">No matched examples</span>'
    parts = ['<div class="examples">']
    for example in examples:
        parts.append('<div class="example">')
        parts.append(f'<code>{html.escape(str(example["labelled_smiles"]))}</code>')
        parts.append(
            '<span>' + html.escape(str(example['identifier'])) +
            f' activity {fmt(example["activity"])} hits {example["nhits"]}</span>'
        )
        parts.append('</div>')
    parts.append('</div>')
    return ''.join(parts)


def build_html(data: dict[str, object], title: str) -> str:
    metadata = data['metadata']
    global_data = data['global']
    rows = data['queries']
    parts = ["""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
"""]
    parts.append(f'<title>{html.escape(title)}</title>\n')
    parts.append("""<style>
:root { --bg: #f6f7f9; --panel: #fff; --ink: #20242a; --muted: #69707d; --border: #d7dce3; --matched: #246b73; --unmatched: #b7c2cf; --warn: #9b5b00; }
* { box-sizing: border-box; }
body { margin: 0; font-family: Arial, Helvetica, sans-serif; background: var(--bg); color: var(--ink); font-size: 14px; }
header { padding: 18px 24px 12px; background: var(--panel); border-bottom: 1px solid var(--border); }
h1 { margin: 0 0 6px; font-size: 22px; line-height: 1.2; }
.summary { margin: 0; color: var(--muted); line-height: 1.4; }
.toolbar { position: sticky; top: 0; z-index: 20; display: flex; flex-wrap: wrap; gap: 10px; align-items: center; padding: 10px 24px; background: #edf1f4; border-bottom: 1px solid var(--border); }
.toolbar input { width: min(420px, 100%); padding: 7px 9px; border: 1px solid var(--border); border-radius: 4px; font: inherit; }
.toolbar button { border: 1px solid var(--border); background: #fff; border-radius: 4px; padding: 7px 10px; font: inherit; cursor: pointer; }
.status { margin-left: auto; color: var(--muted); }
main { padding: 16px 24px 32px; }
.legend { display: flex; gap: 14px; margin-bottom: 12px; color: var(--muted); }
.legend span::before { content: ''; display: inline-block; width: 12px; height: 12px; margin-right: 5px; vertical-align: -1px; }
.legend .m::before { background: var(--matched); }
.legend .u::before { background: var(--unmatched); }
.table-wrap { overflow-x: auto; border: 1px solid var(--border); background: var(--panel); }
table { width: 100%; border-collapse: collapse; min-width: 1180px; }
th, td { border-bottom: 1px solid var(--border); padding: 8px 10px; text-align: left; vertical-align: top; }
th { background: #dfe6eb; color: #3d4651; font-size: 12px; text-transform: uppercase; letter-spacing: 0; white-space: nowrap; }
th.sortable { cursor: pointer; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.name strong { display: block; }
.name code { display: block; color: var(--muted); font-size: 12px; margin-top: 3px; overflow-wrap: anywhere; }
.hist { width: 220px; height: 54px; }
.hist line { stroke: #c7cdd4; stroke-width: 1; }
.hist .matched { fill: var(--matched); }
.hist .unmatched { fill: var(--unmatched); }
.global-hist rect { fill: #7a8794; }
.global-panel { display: inline-flex; align-items: center; gap: 10px; margin-bottom: 12px; color: var(--muted); }
.count-profiles { display: grid; gap: 2px; color: var(--muted); font-size: 12px; margin-top: 4px; }
.warn { color: var(--warn); font-weight: 600; }
.muted { color: var(--muted); }
.examples { display: grid; grid-template-columns: repeat(3, minmax(160px, 1fr)); gap: 6px; max-width: 640px; }
.example { border: 1px solid #e2e5e9; background: #fff; padding: 6px; }
.example code { display: block; overflow-wrap: anywhere; }
.example span { display: block; color: var(--muted); font-size: 12px; margin-top: 4px; }
.hidden { display: none; }
@media (max-width: 720px) { header, .toolbar, main { padding-left: 12px; padding-right: 12px; } .status { width: 100%; margin-left: 0; } }
</style>
</head>
<body>
""")
    parts.append('<header>')
    parts.append(f'<h1>{html.escape(title)}</h1>')
    parts.append(
        '<p class="summary">' +
        f'{metadata["queries"]} queries, {metadata["molecules_read_with_activity"]} molecules with activity. ' +
        f'Activity file: <code>{html.escape(str(metadata["activity"]))}</code>. ' +
        f'Global activity median {fmt(global_data["summary"]["median"])}.</p>'
    )
    parts.append('</header>')
    parts.append("""<div class="toolbar">
<input id="search" type="search" placeholder="Filter by query name or source" oninput="updateView()">
<button type="button" onclick="sortTable('delta_median')">Sort delta median</button>
<button type="button" onclick="sortTable('matched_n')">Sort support</button>
<button type="button" onclick="sortTable('multi')">Sort multiple hits</button>
<span id="status" class="status"></span>
</div>
<main>
""")
    parts.append('<div class="global-panel"><span>Global activity distribution</span>')
    parts.append(simple_histogram_svg(global_data['histogram']))
    parts.append('</div>')
    parts.append('<div class="legend"><span class="m">Matched</span><span class="u">Unmatched</span></div>')
    parts.append("""<div class="table-wrap">
<table id="results">
<thead><tr>
<th>Query</th><th class="sortable" onclick="sortTable('matched_n')">Matched N</th><th>Matched</th><th>Unmatched</th><th class="sortable" onclick="sortTable('delta_median')">Delta Median</th><th class="sortable" onclick="sortTable('delta_mean')">Delta Mean</th><th class="sortable" onclick="sortTable('multi')">Multiple</th><th>Distribution</th><th>Examples</th>
</tr></thead>
<tbody>
""")

    for row in rows:
        matched = row['matched']['summary']
        unmatched = row['unmatched']['summary']
        multi = row['multiplicity']['molecules_with_multiple_matches']
        max_hits = row['multiplicity']['max_matches']
        searchable = f"{row['name']} {row['source']}".lower()
        parts.append(
            f'<tr data-search="{html.escape(searchable, quote=True)}" '
            f'data-matched_n="{matched["n"]}" data-delta_median="{row["delta_median"] if row["delta_median"] is not None else -1e99}" '
            f'data-delta_mean="{row["delta_mean"] if row["delta_mean"] is not None else -1e99}" data-multi="{multi}">'
        )
        parts.append('<td class="name">')
        parts.append(f'<strong>{html.escape(row["name"])}</strong>')
        parts.append(f'<code>{html.escape(row["source"])}</code>')
        parts.append('</td>')
        parts.append(f'<td class="num">{matched["n"]}</td>')
        parts.append(f'<td class="num">median {fmt(matched["median"])}<br>mean {fmt(matched["mean"])}</td>')
        parts.append(f'<td class="num">median {fmt(unmatched["median"])}<br>mean {fmt(unmatched["mean"])}</td>')
        parts.append(f'<td class="num">{fmt(row["delta_median"])}</td>')
        parts.append(f'<td class="num">{fmt(row["delta_mean"])}</td>')
        multi_text = f'{multi}' if multi else '<span class="muted">0</span>'
        if multi:
            multi_text += f'<br><span class="muted">max {max_hits}</span>'
        if metadata.get('multiplicity_mode') == 'count':
            multi_text += multiplicity_profiles_html(row)
        parts.append(f'<td class="num">{multi_text}</td>')
        parts.append('<td>' + histogram_svg(row['matched']['histogram'], row['unmatched']['histogram']) + '</td>')
        parts.append('<td>' + examples_html(row['examples']) + '</td>')
        parts.append('</tr>\n')
    parts.append("""</tbody>
</table>
</div>
</main>
<script>
let sortDescending = true;
function updateView() {
  const needle = document.getElementById('search').value.toLowerCase();
  let shown = 0;
  document.querySelectorAll('#results tbody tr').forEach(row => {
    const match = row.dataset.search.includes(needle);
    row.classList.toggle('hidden', !match);
    if (match) shown += 1;
  });
  document.getElementById('status').textContent = shown + ' queries shown';
}
function sortTable(key) {
  const tbody = document.querySelector('#results tbody');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  rows.sort((a, b) => {
    const av = Number(a.dataset[key]);
    const bv = Number(b.dataset[key]);
    return sortDescending ? bv - av : av - bv;
  });
  sortDescending = !sortDescending;
  rows.forEach(row => tbody.appendChild(row));
  updateView();
}
updateView();
</script>
</body>
</html>
""")
    return ''.join(parts)


def write_json(data: dict[str, object], path: str | None) -> None:
    text = json.dumps(data, indent=2, sort_keys=True) + '\n'
    if path is None or path == '-':
        sys.stdout.write(text)
    else:
        Path(path).write_text(text)


def command_analyse(args: argparse.Namespace) -> int:
    data = analyse(args)
    write_json(data, args.output)
    if args.html:
        Path(args.html).write_text(build_html(data, args.title))
    return 0


def command_html(args: argparse.Namespace) -> int:
    try:
        data = json.loads(Path(args.input).read_text())
    except OSError as exc:
        die(f"Cannot read '{args.input}': {exc}")
    except json.JSONDecodeError as exc:
        die(f"Cannot parse JSON '{args.input}': {exc}")
    Path(args.output).write_text(build_html(data, args.title))
    return 0


def add_analyse_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('--config', required=True, help='SubstructureExploration textproto')
    parser.add_argument('-A', '--activity', required=True, help='Activity file with header, id/activity columns')
    parser.add_argument('-o', '--output', default='-', help='Output JSON file, or - for stdout')
    parser.add_argument('--html', help='Optional HTML report written from the same analysis')
    parser.add_argument('--title', default='Substructure Activity Exploration', help='HTML report title')
    parser.add_argument('--bins', type=int, default=20, help='Number of common histogram bins')
    parser.add_argument('--max_examples', type=int, default=3, help='Matched example molecules retained per query')
    parser.add_argument('--okmissing_activity', action='store_true', help='Allow molecules with no activity value; they are skipped')
    parser.add_argument('--multiplicity', choices=['presence', 'count'], default='presence', help='Primary multiplicity interpretation')
    parser.add_argument('--count-buckets', default='0,1,2,3+', help='Hit count buckets stored in JSON')
    parser.add_argument('smiles', help='Input SMILES file')


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest='command')

    analyse_parser = subparsers.add_parser('analyse', help='Analyse SMILES/activity/query data')
    add_analyse_args(analyse_parser)
    analyse_parser.set_defaults(func=command_analyse)

    html_parser = subparsers.add_parser('html', help='Convert analysis JSON to a static HTML report')
    html_parser.add_argument('input', help='Analysis JSON')
    html_parser.add_argument('-o', '--output', required=True, help='Output HTML file')
    html_parser.add_argument('--title', default='Substructure Activity Exploration', help='HTML report title')
    html_parser.set_defaults(func=command_html)

    args = parser.parse_args(argv)
    if not hasattr(args, 'func'):
        parser.print_help(sys.stderr)
        return 1
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
