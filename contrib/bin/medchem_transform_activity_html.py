#!/usr/bin/env python3
"""Generate an offline HTML report from medchem transformation records.

Run via $LILLYMOL_HOME/run_python.sh. LillyMol handles structural analysis;
RDKit is used only for drawing the selected examples.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from functools import lru_cache
import hashlib
import html
import json
import math
import re
from pathlib import Path
import statistics
import sys

import lillymol

SEPARATOR = '__MCTA__'


def records(path):
    with Path(path).open() as reader:
        for number, line in enumerate(reader, 1):
            if line.strip():
                yield number, line.split()


def number(text, context):
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f'{context}: non-finite activity')
    return value


def product_name(text):
    parts = text.split(SEPARATOR)
    if len(parts) != 2 or not all(parts):
        raise ValueError(f'Invalid product name {text!r}; expected id{SEPARATOR}reaction')
    return tuple(parts)


def read_activity(path):
    if Path(path).suffix.lower() == '.csv':
        with Path(path).open(newline='') as reader:
            rows = [(i, row) for i, row in enumerate(csv.reader(reader), 1) if row]
    else:
        rows = list(records(path))
    if not rows:
        raise ValueError(f'{path}: missing activity header')
    result = {}
    for line, fields in rows[1:]:
        if len(fields) != 2:
            raise ValueError(f'{path}:{line}: expected ID and activity')
        name, value = fields
        if not name or name in result:
            raise ValueError(f'{path}:{line}: empty or duplicate activity ID {name!r}')
        result[name] = number(value, f'{path}:{line}')
    if not result:
        raise ValueError(f'{path}: no activity data')
    return result


def read_starts(path):
    result = {}
    for line, fields in records(path):
        if len(fields) != 2:
            raise ValueError(f'{path}:{line}: expected SMILES and ID')
        smiles, name = fields
        if name in result or SEPARATOR in name:
            raise ValueError(f'{path}:{line}: duplicate ID or reserved separator: {name!r}')
        result[name] = smiles
    if not result:
        raise ValueError(f'{path}: no starting structures')
    return result


@lru_cache(maxsize=4096)
def structure_info(smiles):
    """Only LillyMol determines atom counts and structural identity."""
    molecule = lillymol.MolFromSmiles(smiles)
    if molecule is None or not molecule.natoms():
        raise ValueError(f'LillyMol cannot parse SMILES {smiles!r}')
    return molecule.natoms(), molecule.unique_smiles()


def depiction(smiles):
    """RDKit is a rendering backend, never an analysis/identity backend."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return '<div class="caption">Drawing unavailable<br>' + html.escape(smiles) + '</div>'
    rdDepictor.Compute2DCoords(molecule)
    drawer = rdMolDraw2D.MolDraw2DSVG(230, 165)
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().split('?>', 1)[-1]


def percentile(values, fraction):
    position = (len(values) - 1) * fraction
    index = int(position)
    return values[index] + (values[min(index + 1, len(values) - 1)] - values[index]) * (position - index)


def summarize(values, parents, extent, neutral_width):
    bins = [0] * 40
    if not values:
        return dict(n=0, parents=0, hist=bins, dots=[])
    values.sort()
    for value in values:
        bucket = min(39, max(0, int((value / extent + 1) * 20)))
        bins[bucket] += 1
    return dict(n=len(values), parents=len(parents), mean=statistics.fmean(values),
                median=statistics.median(values), q25=percentile(values, .25),
                q75=percentile(values, .75), min=values[0], max=values[-1],
                positive=sum(v > neutral_width for v in values) / len(values),
                negative=sum(v < -neutral_width for v in values) / len(values),
                hist=bins, dots=values if len(values) <= 150 else [])


def build_data(starting_smiles, activity_file, found=None, notfound=None, predictions=None,
               max_examples=5, neutral_width=0, products_sorted=False):
    starts = read_starts(starting_smiles)
    activity = read_activity(activity_file)
    groups = {}
    prediction_stats = {}
    product_identity = {}
    source_seen = defaultdict(int)
    unscored_products = 0

    def validate(parent):
        if parent not in starts or parent not in activity:
            raise ValueError(f'Missing starting structure or activity for {parent!r}')

    def group(reaction):
        return groups.setdefault(reaction, dict(observed=[], predicted=[], observed_parents=set(),
                                               predicted_parents=set(), examples=[]))

    def example(parent, reaction, smiles, source, target=None, value=None):
        g = group(reaction)
        if products_sorted and source_seen[(source, reaction)] >= max_examples:
            return
        atoms, canonical = structure_info(smiles)
        identity = (parent, canonical, source, target)
        if any(e['_identity'] == identity for e in g['examples']):
            return
        source_seen[(source, reaction)] += 1
        g['examples'].append(dict(parent=parent, before=starts[parent], after=smiles,
                           source=source, product_id=target, atoms=atoms, _identity=identity,
                           starting_activity=activity[parent], product_activity=value,
                           delta=None if value is None else value - activity[parent]))
        # Keep a bounded number from EACH source until the final global merge.
        candidates = sorted([e for e in g['examples'] if e['source'] == source],
                            key=lambda e: (e['atoms'], e['parent'], e['after'], e['product_id'] or ''))
        g['examples'] = [e for e in g['examples'] if e['source'] != source] + candidates[:max_examples]

    if found:
        for line, fields in records(found):
            if len(fields) == 3:
                smiles, combined, target = fields
                parent, reaction = product_name(combined)
            elif len(fields) == 4:
                smiles, parent, reaction, target = fields
            else:
                raise ValueError(f'{found}:{line}: expected SMILES, product name, matched ID')
            validate(parent)
            if target not in activity:
                raise ValueError(f'{found}:{line}: missing activity for {target!r}')
            g = group(reaction)
            delta = activity[target] - activity[parent]
            if not math.isfinite(delta):
                raise ValueError(f'{found}:{line}: non-finite activity change')
            g['observed'].append(delta)
            g['observed_parents'].add(parent)
            example(parent, reaction, smiles, 'observed', target, activity[target])

    if predictions:
        if not notfound:
            raise ValueError('--predictions requires --notfound')
        saw_header = False
        for line, fields in records(predictions):
            if not saw_header:
                saw_header = True
                continue
            if len(fields) == 2:
                parent, reaction = product_name(fields[0])
                value = fields[1]
            elif len(fields) == 3:
                parent, reaction, value = fields
            else:
                raise ValueError(f'{predictions}:{line}: expected product name and prediction')
            validate(parent)
            value = number(value, f'{predictions}:{line}')
            g = group(reaction)
            delta = value - activity[parent]
            if not math.isfinite(delta):
                raise ValueError(f'{predictions}:{line}: non-finite activity change')
            g['predicted'].append(delta)
            g['predicted_parents'].add(parent)
            count, total = prediction_stats.get((parent, reaction), (0, 0.0))
            prediction_stats[(parent, reaction)] = count + 1, total + value
        if not saw_header:
            raise ValueError(f'{predictions}: missing header')
        for line, fields in records(notfound):
            if len(fields) == 2:
                smiles, combined = fields
                parent, reaction = product_name(combined)
            elif len(fields) == 3:
                smiles, parent, reaction = fields
            else:
                raise ValueError(f'{notfound}:{line}: expected SMILES and product name')
            validate(parent)
            key = parent, reaction
            if key not in prediction_stats:
                unscored_products += 1
                continue
            first, ambiguous = product_identity.get(key, (smiles, False))
            if first != smiles and not ambiguous:
                ambiguous = structure_info(first)[1] != structure_info(smiles)[1]
            product_identity[key] = first, ambiguous
            example(parent, reaction, smiles, 'predicted')
        missing = prediction_stats.keys() - product_identity.keys()
        if missing:
            raise ValueError(f'Predictions have no notfound structure: {sorted(missing)[:3]}')

    if not groups:
        raise ValueError('No transformation evidence found')
    extent = max(abs(v) for g in groups.values() for source in ['observed', 'predicted'] for v in g[source]) or 1.0
    drawings = {}

    def image(smiles):
        structure_info(smiles)  # Validate selected starts with LillyMol too.
        key = hashlib.sha256(smiles.encode()).hexdigest()
        if key not in drawings:
            drawings[key] = depiction(smiles)
        return key

    rows = []
    for reaction, g in sorted(groups.items()):
        examples = sorted(g['examples'], key=lambda e: (e['atoms'], e['parent'], e['after'], e['source'], e['product_id'] or ''))[:max_examples]
        for e in examples:
            if e['source'] == 'predicted':
                key = e['parent'], reaction
                e['ambiguous'] = product_identity[key][1]
                count, total = prediction_stats[key]
                e['prediction_count'] = count
                if not e['ambiguous']:
                    e['product_activity'] = total / count
                    e['delta'] = total / count - e['starting_activity']
            e.pop('_identity')
            e['before_svg'] = image(e['before'])
            e['after_svg'] = image(e['after'])
        rows.append(dict(name=reaction, examples=examples,
                         observed=summarize(g['observed'], g['observed_parents'], extent, neutral_width),
                         predicted=summarize(g['predicted'], g['predicted_parents'], extent, neutral_width)))
    return dict(rows=rows, drawings=drawings, extent=extent, neutral_width=neutral_width,
                max_examples=max_examples, unscored_products=unscored_products)


def render_html(data, title):
    # Escape '<' so names cannot terminate the inert JSON script element.
    encoded = json.dumps(data, allow_nan=False).replace('<', '\\u003c').replace('\u2028', '\\u2028').replace('\u2029', '\\u2029')
    template = Path(__file__).with_suffix('.template.html').read_text()
    replacements = {'TITLE': html.escape(title), 'DATA': encoded}
    return re.sub(r'@@(TITLE|DATA)@@', lambda match: replacements[match.group(1)], template)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--starting-smiles', required=True, type=Path)
    parser.add_argument('--activity', '-A', required=True, type=Path)
    parser.add_argument('--found', type=Path)
    parser.add_argument('--notfound', type=Path)
    parser.add_argument('--predictions', type=Path)
    parser.add_argument('--max-examples', type=int, default=5)
    parser.add_argument('--neutral-width', type=float, default=0)
    parser.add_argument('--products-sorted', action='store_true', help='Products already sorted by increasing LillyMol atom count')
    parser.add_argument('--title', default='Medchem transformation activity')
    parser.add_argument('--output', '-o', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.max_examples < 1 or not math.isfinite(args.neutral_width) or args.neutral_width < 0:
        parser.error('--max-examples must be positive; --neutral-width must be finite and nonnegative')
    if not args.found and not args.predictions:
        parser.error('Specify --found and/or --predictions')
    if args.notfound and not args.predictions:
        parser.error('--notfound requires --predictions')
    try:
        for path in [args.starting_smiles, args.activity, args.found, args.notfound, args.predictions]:
            if path and (args.output.resolve() == path.resolve() or
                         (args.output.exists() and args.output.samefile(path))):
                raise ValueError(f'Output would overwrite input {path}')
        data = build_data(args.starting_smiles, args.activity, args.found, args.notfound,
                          args.predictions, args.max_examples, args.neutral_width, args.products_sorted)
        args.output.write_text(render_html(data, args.title), encoding='utf-8')
        print(f'Wrote {len(data["rows"])} transformations to {args.output}; '
              f'{data["unscored_products"]} notfound products had no prediction', file=sys.stderr)
    except (OSError, ValueError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
