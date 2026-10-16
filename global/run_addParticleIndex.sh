#!/usr/bin/env bash
set -euo pipefail # Exit on error, undefined variable, and fail on pipe errors

if [ $# -ne 4 ] || ! [[ $1 =~ ^(0|[1-9][0-9]*)$ ]] || ! [[ $2 =~ ^(0|[1-9][0-9]*)$ ]]; then  # check if integers
    echo "Usage: $0 <N> <M> <onSTBC> <path>  (run add_particle_index_to_tracks.C for <path>/seed_N..M)"
    exit 1
fi

N=$1
M=$2
onSTBC=$3
path=$4

for n in $(seq "$N" "$M"); do
    full_path="${path}/seed_${n}"
    echo "=== Processing ${full_path} ==="
    root -l -b -q "add_particle_index_to_tracks.C(\"${full_path}\", $onSTBC)"
done

echo "Done."
