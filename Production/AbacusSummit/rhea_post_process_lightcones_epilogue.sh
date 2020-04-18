#!/usr/bin/env bash

# Called with the lightcones.concat dir to process.
# For each subdir, do:
# - Merge checksum files
# - Copy in headers
# - Zip headers
# Then, delete original lightcones. Rename new lightcones.

set -e

NEWLCDIR=$1
echo "Starting lightcones epilogue in $NEWLCDIR"
cd $NEWLCDIR

LCDIR=$(realpath ../lightcone)

for D in heal/ pid/ rv/; do
    pushd $D

    for S in $LCDIR/Step*; do
        cp $S/header header_$(basename $S)
    done

    zip --move -qr headers.zip header_*
    $ABACUS/external/fast-cksum/bin/fast_cksum headers.zip > headers.zip.crc32

    $ABACUS/external/fast-cksum/bin/merge_checksum_files.py --delete *.crc32 > checksums.crc32

    popd > /dev/null
done

# Back to SIMDIR
cd ..

echo "Deleting lightcone/"
time rm -rf lightcone/

echo "Renaming lightcones.concat/"
mv lightcones.concat/ lightcones/
