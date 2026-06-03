#!/usr/bin/env bash
set -euo pipefail # Exit on error, undefined variable, and fail on pipe errors

if [ $# -ne 1 ] || ! [[ $1 =~ ^[1-9][0-9]*$ ]]; then  # check if integer
    echo "Usage: $0 <N>   (run addMatchingToParticles.C for pythia_20k_1 .. pythia_20k_N)"
    exit 1
fi

N=$1

for n in $(seq 1 "$N"); do
    path="pythia_20k/pythia_20k_${n}"
    echo "=== Processing ${path} ==="
    root -l -b -q "addMatchingToParticles.C(\"${path}\")"
done

echo "Done."
