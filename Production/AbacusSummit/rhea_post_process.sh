#!/usr/bin/env bash

set -e  # exit on error
set -o pipefail

export SIM_SET=${SIM_SET:-AbacusSummit}
SIM_NAME=$1

LOGDIR=$(pwd)/logs/$SIM_NAME
mkdir -p $LOGDIR

##### queue up group processing #####

mkdir -p tmp
DISBATCH_TASKFILE=$(pwd)/tmp/${SIM_NAME}.group.disbatch

CHUNK=50
NWORKER=8
GROUPSCRIPT="$ABACUS/Abacus/convert_raw_groups_to_asdf.py --chunk=$CHUNK --nworkers=$NWORKER"
GROUPDIR=$ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME/group

echo "#DISBATCH PREFIX cd $GROUPDIR; $GROUPSCRIPT " > $DISBATCH_TASKFILE  # write
echo "#DISBATCH SUFFIX > $(pwd)/\${DISBATCH_NAMETASKS}_\${DISBATCH_JOBID}_\${DISBATCH_TASKID}.log 2>&1" >> $DISBATCH_TASKFILE  # append

(cd $GROUPDIR; echo Step*/) | tr " " "\n" >> $DISBATCH_TASKFILE

mkdir -p ./logs/$SIM_NAME/disbatchGroups

# We can fill up to ~36 nodes
NNODES=3
WALLTIME=3:00:00
# 1 task per node. Each task uses 8 workers, each running 2 blosc threads
JOBID1=$(sbatch -t $WALLTIME -N $NNODES --ntasks-per-node=1 \
        --mem=0 -p batch -A AST145 --parsable --job-name=${SIM_NAME}_PostProcessGroups -o logs/$SIM_NAME/%x.out \
        --wrap "./disBatch/disBatch.py -e -p ./logs/$SIM_NAME/disbatchGroups/$SIM_NAME $DISBATCH_TASKFILE && echo Groups completed successfully.")
EPILOG_DEPEND+=":$JOBID1"

##### queue up lightcone processing #####

LCDIR=$ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME/lightcone/

if [ -d $LCDIR ]; then

    # Write a disBatch file with one line per light cone per file type per step
    DISBATCH_TASKFILE=$(pwd)/tmp/${SIM_NAME}.LC.disbatch

    NEWLCDIR=$ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME/lightcones.concat/

    LCSCRIPT=$(pwd)/rhea_post_process_lightcones.py
    echo "#DISBATCH PREFIX cd $LCDIR; $LCSCRIPT " > $DISBATCH_TASKFILE  # write
    echo "#DISBATCH SUFFIX > $(pwd)/\${DISBATCH_NAMETASKS}_\${DISBATCH_JOBID}_\${DISBATCH_TASKID}.log 2>&1 " >> $DISBATCH_TASKFILE  # append

    # Args could get dangerously long, use relative paths
    (cd $LCDIR; find Step*/ -name 'LightCone*' | cut -d. -f1 | sort -u) >> $DISBATCH_TASKFILE  # Step0500/LightCone0_heal, etc

    # Now generate the checksum concatenation tasks

    LCEPI_SCRIPT=$(pwd)/rhea_post_process_lightcones_epilogue.sh
    cat >> $DISBATCH_TASKFILE << EOM
#DISBATCH BARRIER CHECK
#DISBATCH PREFIX cd $NEWLCDIR/
#DISBATCH SUFFIX ; $LCEPI_SCRIPT > $(pwd)/\${DISBATCH_NAMETASKS}_\${DISBATCH_JOBID}_\${DISBATCH_TASKID}.log 2>&1 
heal/
pid/
rv/
EOM

    cat >> $DISBATCH_TASKFILE << EOM
#DISBATCH BARRIER CHECK
#DISBATCH PREFIX 
#DISBATCH SUFFIX 
#rm -rf $LCDIR
##DISBATCH BARRIER CHECK
mv $NEWLCDIR $(dirname $LCDIR)/lightcones
EOM

    # We can fill up to ~2000 nodes
    NNODES=3
    WALLTIME=1:00:00

    mkdir -p ./logs/$SIM_NAME/disbatchLC

    # 4 cores per task: only two blosc threads per task, but cores have two hyperthreads.  So 8 tasks per node.
    JOBID2=$(sbatch -t $WALLTIME -N $NNODES -c 4 \
            -o logs/$SIM_NAME/%x.out --mem=0 -A AST145 -p batch --parsable --job-name=${SIM_NAME}_PostProcessLC \
            --wrap "./disBatch/disBatch.py -e -p ./logs/$SIM_NAME/disbatchLC/$SIM_NAME $DISBATCH_TASKFILE && echo Lightcones completed successfully.")
    EPILOG_DEPEND+=":$JOBID2"
fi

##### Queue up the time slice processing #####

SIMDIR=$ABACUSSUMMIT_PERSIST/$SIM_SET/$SIM_NAME
SLICEDIRS=$(cd $SIMDIR; find -maxdepth 1 -name 'slice*.*')

if [[ -n "$SLICEDIRS" ]]; then

    # Write a disBatch file with one line per chunk per file type per time slice
    # Chunks should be specified with a literal pattern like slice1.000/*slab001?.L0_pack9.dat
    DISBATCH_TASKFILE=$(pwd)/tmp/${SIM_NAME}.timeslice.disbatch

    TSSCRIPT=$(pwd)/rhea_post_process_timeslice.py
    echo "#DISBATCH PREFIX cd ${SIMDIR}; ${TSSCRIPT} \"" > $DISBATCH_TASKFILE  # write
    echo "#DISBATCH SUFFIX \"" >> $DISBATCH_TASKFILE  # append

    for SLICE in $SLICEDIRS; do
        # Write one line for every chunk of 10 for every file type
        # slice1.000/*slab001?.L0_pack9.dat

        i=0
        while : ; do
            FNPREFIX="${SLICE}/$(printf '*slab%03d?*.dat' $i)"
            MATCHES=$(cd $SIMDIR; find -wholename "$FNPREFIX" -print -quit)
            if [[ -z "$MATCHES" ]]; then
                break
            fi
            for FTYPE in L0_pack9 L0_pack9_pids field_pack9 field_pack9_pids; do
                echo "${SLICE}"/$(printf "*slab%03d?.${FTYPE}.dat" $i) >> $DISBATCH_TASKFILE
            done
            i=$((i + 1))
        done
    done

    # Now generate the checksum concatenation tasks
    TSEPI_SCRIPT=$(pwd)/rhea_post_process_timeslice_epilogue.sh
    cat >> $DISBATCH_TASKFILE << EOM
#DISBATCH BARRIER CHECK
#DISBATCH PREFIX 
#DISBATCH SUFFIX 
$TSEPI_SCRIPT $SIMDIR/slices/z*
EOM

#    cat >> $DISBATCH_TASKFILE << EOM
##DISBATCH BARRIER CHECK
#rm -rf $SIMDIR/slice*.*
#EOM

    # We can fill up to ??? nodes
    NNODES=3
    WALLTIME=1:00:00

    mkdir -p ./logs/$SIM_NAME/disbatchTS

    JOBID3=$(sbatch -t $WALLTIME -N $NNODES -c 4 \
            -o logs/$SIM_NAME/%x.out --mem=0 -A AST145 -p batch --parsable --job-name=${SIM_NAME}_PostProcessTS \
            --wrap "./disBatch/disBatch.py -e -p ./logs/$SIM_NAME/disbatchTS/$SIM_NAME $DISBATCH_TASKFILE && echo Time slices completed successfully.")
    EPILOG_DEPEND+=":$JOBID3"
fi

##### Queue up the epilogue #####
# This job will clean up and mark completion after all previous jobs have completed
# It carries the special job name that the assembly line will look for to determine if post-processing is queued
sbatch --job-name=${SIM_NAME}_PostProcessEpilogue -o logs/$SIM_NAME/%x.out --depend=afterok${EPILOG_DEPEND} --kill-on-invalid-dep=yes rhea_post_process_epilogue.slurm $SIM_NAME
