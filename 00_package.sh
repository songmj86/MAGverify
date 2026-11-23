#!/usr/bin/env bash
# Create an archive whose contents extract directly into the script directory.
set -euo pipefail

repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
archive="$(dirname -- "$repo_dir")/script_github.tar.gz"

tar -czf "$archive" \
  --exclude='.git' \
  --exclude='__pycache__' \
  --exclude='*.pyc' \
  --exclude='.pytest_cache' \
  -C "$repo_dir" .

printf 'Created: %s\n' "$archive"
