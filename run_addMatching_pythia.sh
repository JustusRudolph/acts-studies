#!/usr/bin/env bash
set -euo pipefail # Exit on error, undefined variable, and fail on pipe errors

if [ $# -ne 1 ] || ! [[ $1 =~ ^[1-9][0-9]*$ ]]; then  # check if integer
    echo "Usage: $0 <N>   (run addMatchingToParticles.C for 50k_pythia_1 .. 50k_pythia_N)"
    exit 1
fi

N=$1

for n in $(seq 1 "$N"); do
    path="50k_pythia_${n}"
    echo "=== Processing ${path} ==="
    root -l -b -q "addMatchingToParticles.C(\"${path}\")"
done

echo "Done."
