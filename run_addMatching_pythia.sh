#!/usr/bin/env bash
set -euo pipefail # Exit on error, undefined variable, and fail on pipe errors

if [ $# -ne 3 ] || ! [[ $1 =~ ^[1-9][0-9]*$ ]] || ! [[ $2 =~ ^[1-9][0-9]*$ ]]; then  # check if integers
    echo "Usage: $0 <N> <M> <onSTBC>  (run add_matching_to_particles.C for pythia_2500ev/seed_N..M)"
    exit 1
fi

N=$1
M=$2
onSTBC=$3

for n in $(seq "$N" "$M"); do
    path="pythia_2500ev/seed_${n}"
    echo "=== Processing ${path} ==="
    root -l -b -q "add_matching_to_particles.C(\"${path}\", $onSTBC)"
done

echo "Done."
