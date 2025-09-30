#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   sh ./run_chunks.sh /path/to/inputs.list your_fcl.fcl NUM_CHUNKS /path/to/outdir
#
# Each chunk runs 10,000 events, increases --nskip by 10,000, and
# writes numbered event (-o) and histogram (-T) outputs.

FILELIST="${1:?need file list}"
FCL="${2:?need .fcl file}"
NCHUNKS="${3:?need number of chunks}"
OUTDIR="${4:?need output directory}"
OUTNAME="${5:?need output name prefix}"

# Basic checks
[[ -f "$FILELIST" ]] || { echo "File list not found: $FILELIST" >&2; exit 1; }
[[ -s "$FILELIST" ]] || { echo "File list is empty: $FILELIST" >&2; exit 1; }

mkdir -p "$OUTDIR"

for (( i=0; i< NCHUNKS; ++i )); do
  idx=$(printf "%05d" "$i")
  N=10000
  NSKIP=$(( i * N ))

  HISTO_OUT="${OUTDIR}/${OUTNAME}_histo_${idx}.root"
  EVENT_OUT="${OUTDIR}/${OUTNAME}_out_${idx}.root"

  echo "=== Chunk $i ==="
  echo "  -S $FILELIST"
  echo "  -n $N --nskip $NSKIP"
  echo "  -T $HISTO_OUT"
  echo "  -o $EVENT_OUT"

  lar -c "$FCL" \
      -S "$FILELIST" \
      -n "$N" \
      --nskip "$NSKIP" \
      -T "$HISTO_OUT" \
      -o "$EVENT_OUT"
done
