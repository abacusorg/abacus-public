#!/usr/bin/env bash

# Called with each new timeslice dir. Do:
# - Copy in headers
# - Merge checksum files

set -e

ALL_ZDIRS=$@

for ZDIR in $ALL_ZDIRS; do 
    SIM_DIR=$ZDIR/../..
    SLICEDIR=$(basename $ZDIR)
    SLICEDIR=$SIM_DIR/${SLICEDIR/z/slice}

    # Copy header to redshift dir and each type subdir for completeness
    cp $SLICEDIR/header $ZDIR
    for TYPEDIR in $ZDIR/*/; do
        cd $TYPEDIR
        cp $SLICEDIR/header .
        $ABACUS/external/fast-cksum/bin/fast_cksum header > header.crc32

        $ABACUS/external/fast-cksum/bin/merge_checksum_files.py --delete *.crc32 > checksums.crc32
        cd - > /dev/null
    done
done
