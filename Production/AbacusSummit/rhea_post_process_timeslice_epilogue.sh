#!/usr/bin/env bash

# Called with new slices dir.  For each z subdir, do:
# - Copy in headers
# - Merge checksum files
# Then, delete original slices.

set -e

NEWSLICES=$1
echo "Starting time slice epilogue in $NEWSLICES"
cd $NEWSLICES

SIM_DIR=$(realpath ..)

for ZDIR in z*/; do 
    SLICEDIR=$(basename $ZDIR)
    SLICEDIR=$SIM_DIR/${SLICEDIR/z/slice}

    # Copy header to redshift dir and each type subdir for completeness
    cp $SLICEDIR/header $ZDIR
    for TYPEDIR in $ZDIR/*/; do
        pushd $TYPEDIR
        cp $SLICEDIR/header .
        $ABACUS/external/fast-cksum/bin/fast_cksum header > header.crc32

        $ABACUS/external/fast-cksum/bin/merge_checksum_files.py --delete *.crc32 > checksums.crc32
        popd > /dev/null
    done
done

echo "Deleting slice*.*/"  # careful not to glob slices/!
time rm -rf ${SIM_DIR}/slice*.*/
