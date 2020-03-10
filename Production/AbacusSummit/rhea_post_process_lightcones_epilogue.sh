#!/usr/bin/env bash

# Called from inside each new lightcone dir. Do:
# - Merge checksum files
# - Copy in headers
# - Zip headers

set -e

LCDIR=../../lightcone

for S in $LCDIR/Step*; do
    cp $S/header header_$(basename $S)
done

zip --move -qr headers.zip header_*
$ABACUS/external/fast-cksum/bin/fast_cksum headers.zip > headers.zip.crc32

$ABACUS/external/fast-cksum/bin/merge_checksum_files.py --delete *.crc32 > checksums.crc32
