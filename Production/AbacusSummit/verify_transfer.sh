#!/usr/bin/env bash

# File purges at OLCF may mean that transfers could be missing files.
# But the purges don't delete directories, so we can check that all
# directories (1) have a checksums.crc32 file, and (2) have all files
# listed in that file.

# Usage:
# verify_transfer.sh <SIM_DIR>

# allow **
shopt -s globstar

for DIR in $1/**/; do
    if [[ ! -f $DIR/checksums.crc32 ]]; then
        # These are the directories that are allowed to not have checksums.crc32
        if echo $(basename $DIR) | \grep -Pq '^(halos|z.*|lightcones|slices)$'; then
            continue
        fi
        echo "$DIR missing checksums.crc32"
        continue
    fi
    FNS=$(awk '{print $3}' $DIR/checksums.crc32)
    pushd $DIR
    for FN in $FNS; do
        if [[ ! -f $FN ]]; then
            echo "$FN missing"
        fi
    done
    popd > /dev/null
done
