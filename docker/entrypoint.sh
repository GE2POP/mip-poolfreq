#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -eq 0 ] || [[ "${1:-}" =~ ^(-h|--help|help)$ ]]; then
  cat <<'EOF'
mip-poolfreq container

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
