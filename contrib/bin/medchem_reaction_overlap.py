#!/usr/bin/env python3
"""Compare a proposed MedchemWizard reaction with a library on a SMILES corpus.

Run via ${LILLYMOL_HOME}/run_python.sh; requires the nanobind lillymol module.
All equivalence statements apply only to the supplied input records.
"""

import argparse
import json
import os
from pathlib import Path
import sys

import lillymol


def read_molecules(paths):
    """Use record positions as keys: names need not be unique or single tokens."""
    records = []
    for path in paths:
        with Path(path).open() as inp:
            for line_number, line in enumerate(inp, 1):
                if not line.strip() or line.lstrip().startswith('#'):
                    continue
                mol = lillymol.MolFromSmiles(line.strip())
                if mol is None:
                    raise ValueError(f'{path}:{line_number}: invalid SMILES')
                records.append((str(path), line_number, mol))
    if not records:
        raise ValueError('No input molecules')
    return records


def read_library(path):
    """Preserve the format marker; .rxn alone does not identify the format."""
    entries = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        proto = line.startswith('PROTO:')
        filename = line[6:] if proto else line
        entries.append(((path.parent / filename).resolve(), proto))
    if not entries:
        raise ValueError(f'Empty reaction list: {path}')
    return entries


def load_reaction(path, proto):
    wizard = lillymol.MedchemWizard()
    wizard.read_reaction(str(path), textproto=proto)
    return wizard


def generate_products(wizard, records):
    products = set()
    for index, (_, _, mol) in enumerate(records):
        for product in wizard.process(mol):
            products.add((index, product.unique_smiles()))
    return products


def examples(keys, records, limit):
    result = []
    for index, product in sorted(keys)[:limit]:
        source, line, mol = records[index]
        result.append(dict(source=source, line=line, name=mol.name(),
                           parent=mol.smiles(), product=product))
    return result


def compare(candidate, existing, records, limit):
    shared = candidate & existing
    candidate_only = candidate - existing
    existing_only = existing - candidate
    if not shared:
        relation = 'no_overlap'
    elif not candidate_only and not existing_only:
        relation = 'identical_on_corpus'
    elif not candidate_only:
        relation = 'candidate_contained_on_corpus'
    elif not existing_only:
        relation = 'existing_contained_on_corpus'
    else:
        relation = 'partial_overlap'
    return dict(relation=relation, existing_products=len(existing), shared=len(shared),
                candidate_only=len(candidate_only), existing_only=len(existing_only),
                candidate_overlap=len(shared) / len(candidate) if candidate else 0,
                existing_overlap=len(shared) / len(existing) if existing else 0,
                examples={
                    'shared': examples(shared, records, limit),
                    'candidate_only': examples(candidate_only, records, limit),
                    'existing_only': examples(existing_only, records, limit)})


def check(candidate_path, library, inputs, candidate_format='auto', max_examples=5):
    candidate_path = candidate_path.resolve()
    entries = read_library(library.resolve())
    if candidate_format == 'auto':
        formats = {proto for path, proto in entries if path == candidate_path}
        if len(formats) > 1:
            raise ValueError('Candidate has conflicting formats in the library')
        if formats:
            proto = formats.pop()
        else:
            # New MedchemWizard reactions are normally textproto, even as .rxn.
            proto = True
    else:
        proto = candidate_format == 'proto'
    records = read_molecules(inputs)
    candidate = generate_products(load_reaction(candidate_path, proto), records)
    if not candidate:
        raise ValueError('Proposed reaction produced no products; supply matching molecules')
    rows = []
    errors = []
    seen = set()
    for path, entry_proto in entries:
        if path == candidate_path or (path, entry_proto) in seen:
            continue
        seen.add((path, entry_proto))
        try:
            wizard = load_reaction(path, entry_proto)
        except (OSError, ValueError, RuntimeError) as exc:
            errors.append(dict(reaction=str(path), error=str(exc)))
            continue
        existing = generate_products(wizard, records)
        row = compare(candidate, existing, records, max_examples)
        row.update(reaction=str(path), name=path.name)
        rows.append(row)
    rows.sort(key=lambda row: (-row['shared'], row['candidate_only'] +
                              row['existing_only'], row['reaction']))
    return dict(candidate=str(candidate_path), library=str(library.resolve()),
                input_records=len(records), candidate_products=len(candidate),
                reactions_checked=len(rows), results=rows, errors=errors,
                interpretation='Results describe only the supplied corpus, not universal equivalence.')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reaction', type=Path, help='Proposed reaction file')
    parser.add_argument('smiles', nargs='+', type=Path, help='Representative SMILES files')
    home = Path(os.environ.get('LILLYMOL_HOME', Path(__file__).resolve().parents[2]))
    parser.add_argument('--library', type=Path, default=home / 'data/MedchemWizard/REACTIONS')
    parser.add_argument('--format', choices=('auto', 'proto', 'msi'), default='auto',
                        help='Candidate format: auto uses its library entry, otherwise proto')
    parser.add_argument('--json', type=Path, help='Detailed report with bounded examples')
    parser.add_argument('--max-examples', type=int, default=5,
                        help='Examples per category per reaction (default: 5)')
    args = parser.parse_args(argv)
    if args.max_examples < 0:
        parser.error('--max-examples must be nonnegative')
    try:
        return run(args, parser)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        return 1


def run(args, parser):
    # Never allow an output path to overwrite any source.
    if args.json:
        protected = {p.resolve() for p in [args.reaction, args.library, *args.smiles]}
        protected.update(path for path, _ in read_library(args.library.resolve()))
        if args.json.resolve() in protected:
            parser.error('--json would overwrite an input or reaction file')
    report = check(args.reaction, args.library, args.smiles, args.format,
                   args.max_examples)
    if args.json:
        args.json.write_text(json.dumps(report, indent=2) + '\n')
    print(f"{report['candidate_products']} candidate products from "
          f"{report['input_records']} records; checked {report['reactions_checked']} reactions.")
    print('reaction\trelationship\tshared\tcandidate_only\texisting_only\tcandidate_overlap\texisting_overlap')
    for row in report['results']:
        if row['shared']:
            print(f"{row['name']}\t{row['relation']}\t{row['shared']}\t"
                  f"{row['candidate_only']}\t{row['existing_only']}\t"
                  f"{row['candidate_overlap']:.1%}\t{row['existing_overlap']:.1%}")
    print(report['interpretation'])
    for error in report['errors']:
        print(f"Incomplete comparison: {error['error']}", file=sys.stderr)
    return 2 if report['errors'] else 0


if __name__ == '__main__':
    sys.exit(main())
