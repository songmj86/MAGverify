#!/usr/bin/env python3
"""Stream read files to stdout and report progress with the remaining time on stderr.

Progress is the share of input bytes passed on. The rate is measured from the moment
the consumer (for example minimap2 after building its index) starts reading, so index
building does not distort the remaining-time estimate.
"""
import argparse
import os
import sys
import time

CHUNK = 4 * 1024 * 1024
WARMUP_BYTES = 16 * 1024 * 1024


def duration(seconds):
    seconds = int(max(0, seconds))
    hours, rest = divmod(seconds, 3600)
    minutes = rest // 60
    if hours:
        return f'{hours}h{minutes:02d}m'
    if minutes:
        return f'{minutes}m'
    return f'{seconds}s'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='+')
    parser.add_argument('--label', default='reads')
    parser.add_argument('--interval', type=float, default=300.0, help='Seconds between reports')
    args = parser.parse_args()

    gzipped = set()
    for path in args.files:
        with open(path, 'rb') as handle:
            gzipped.add(handle.read(2) == b'\x1f\x8b')
    if len(gzipped) > 1:
        sys.exit('ERROR: read files mix gzip-compressed and uncompressed files')
    total = sum(os.path.getsize(path) for path in args.files) or 1
    start = time.monotonic()
    last_report = start
    rate_start = None
    rate_bytes = 0
    done = 0
    out = sys.stdout.buffer
    try:
        for path in args.files:
            with open(path, 'rb') as handle:
                while True:
                    block = handle.read(CHUNK)
                    if not block:
                        break
                    out.write(block)
                    done += len(block)
                    now = time.monotonic()
                    if rate_start is None and done >= WARMUP_BYTES:
                        rate_start, rate_bytes = now, done
                    if now - last_report >= args.interval:
                        last_report = now
                        remaining = 'estimating'
                        if rate_start is not None and now > rate_start and done > rate_bytes:
                            speed = (done - rate_bytes) / (now - rate_start)
                            remaining = duration((total - done) / speed)
                        print(f'[{args.label}] {100 * done / total:.0f}% | elapsed '
                              f'{duration(now - start)} | remaining {remaining}',
                              file=sys.stderr, flush=True)
        out.flush()
    except BrokenPipeError:
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        sys.exit(1)
    print(f'[{args.label}] 100% read | elapsed {duration(time.monotonic() - start)}; '
          'finishing alignment and sorting', file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
