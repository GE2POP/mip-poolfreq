#!/usr/bin/env bash
set -euo pipefail

VERSION_FILE="/opt/mip-poolfreq/VERSION"
PROJECT_VERSION="unknown"

if [ -f "$VERSION_FILE" ]; then
  PROJECT_VERSION="$(cat "$VERSION_FILE")"
fi

if [ "$#" -eq 0 ] || [[ "${1:-}" =~ ^(-h|--help|help)$ ]]; then
  cat <<EOF
mip-poolfreq container
Version: ${PROJECT_VERSION}

This container provides two command-line tools:

1. Afreq
   Estimate allele frequencies from sequencing data.

2. Gfreq
   Estimate genotype frequencies in mixtures of known components.

Usage:
  docker run --rm <image> Afreq --help
  docker run --rm <image> Gfreq --help

----- Afreq help -----
EOF
  Afreq --help || true
  echo
  echo "----- Gfreq help -----"
  Gfreq --help || true
  exit 0
fi

exec "$@"
