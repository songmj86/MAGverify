#!/usr/bin/env python3
import argparse
import csv
import os
import tempfile


def parse_csv(text):
    return [x.strip() for x in (text or '').split(',') if x.strip()]


def main():
    p = argparse.ArgumentParser(description='Export key columns plus evidence-specific columns from an integrated TSV.')
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--keys', required=True, help='Comma-separated key/context columns that must exist and are retained.')
    p.add_argument('--prefix', action='append', default=[], help='Column prefix to retain. Repeatable.')
    p.add_argument('--exclude-prefix', action='append', default=[], help='Prefix to exclude after inclusion. Repeatable.')
    p.add_argument('--include', default='', help='Additional comma-separated exact columns to retain when present.')
    args = p.parse_args()

    with open(args.input, encoding='utf-8-sig', newline='') as h:
        r = csv.DictReader(h, delimiter='\t')
        fields = r.fieldnames or []
        rows = list(r)
    if not fields:
        raise SystemExit(f'ERROR: TSV has no header: {args.input}')

    requested_keys = parse_csv(args.keys)
    missing_keys = [x for x in requested_keys if x not in fields]
    if missing_keys:
        raise SystemExit(
            f'ERROR: requested validation key column(s) are absent from {args.input}: ' + ', '.join(missing_keys)
        )
    keys = requested_keys
    exact = [x for x in parse_csv(args.include) if x in fields]
    selected = []
    for col in fields:
        if col in keys or col in exact:
            selected.append(col)
            continue
        if args.prefix and not any(col.startswith(pref) for pref in args.prefix):
            continue
        if args.prefix and any(col.startswith(pref) for pref in args.exclude_prefix):
            continue
        if args.prefix:
            selected.append(col)

    # Preserve order, avoid duplicates.
    out_fields = []
    seen = set()
    for col in selected:
        if col not in seen:
            seen.add(col)
            out_fields.append(col)
    if not out_fields:
        raise SystemExit(f'ERROR: no requested validation columns found in: {args.input}')

    dest = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(dest)+'.', suffix='.tmp', dir=os.path.dirname(dest), text=True)
    os.close(fd)
    try:
        with open(tmp, 'w', encoding='utf-8', newline='') as h:
            w = csv.DictWriter(h, fieldnames=out_fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


if __name__ == '__main__':
    main()
