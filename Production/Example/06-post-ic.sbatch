#!/bin/bash

# Create density and displacement IC files at the given resolution(s)
# Depends on: none

#SBATCH -N 1
#SBATCH -t 0-0:30
#SBATCH -C cpu
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -J post-ic
#SBATCH --mem=0

set -e

simpath=$1
allres="576 1152"

#########

if [ -z "$simpath" ]; then
    echo "Usage: sbatch $0 simpath"
    exit 1
fi

simname=$(basename $simpath)
out=${SCRATCH}/abacus_ic/

mkdir -p ${out}

for res in $allres; do
    scratch=/dev/shm/lgarriso/ic-$simname-$res/
    mkdir -p ${scratch}

    # TODO: parallel safe?
    python -m Abacus.zeldovich --density --ZD_qPLT=0 --ZD_NumBlock=32 \
        --ICFormat=\"ZelSimple\" --CPD=1 --out-parent=${scratch} \
        --ppd=${res} --ZD_CornerModes=1 \
        $simpath/abacus.par
        
    $ABACUS/Production/AbacusSummit/nersc_package_ics.py \
        ${scratch}/${simname}/ic/ --delete --out=${out}
done

# Only do on the last round
$ABACUS/external/fast-cksum/bin/merge_checksum_files.py \
    --delete ${out}/${simname}/*.crc32 \
    | sponge ${out}/${simname}/checksums.crc32
