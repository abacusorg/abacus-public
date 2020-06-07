#!/bin/bash
set -e

# This script runs on the summit compute node
# Fail if any command returns a non-zero exit status

#Specify simulation to run on this node.
######################

export SIM_NAME="AbacusSummit_small_test_${JSM_NAMESPACE_RANK}"
export MINISUITE="MiniTestSuite"
######################

# Specify the directory holding the simulation config inside the abacus/Production dir
######################
SET_NAME="AbacusSummit"
if [[ -z "$SIM_NAME" ]]; then
    echo "Need $SIM_NAME!" ; exit 1
fi
######################

echo -e "* Loading modules\n"  
source modules.sh
echo -e "\n\n\n\n"  

export ABACUS_PERSIST=$BBPATH
export ABACUS_DERIVS=${PROJWORK} # extra env. var. to search for Derivatives on GPFS and copy to BB.
echo -e "ABACUS_DERIVS set to ${ABACUS_DERIVS}"  

######################
###################### # performed automatically by serial abacus. No need to separate out as in rhea/summit production line. 
# Run ICs # 
#PARFN=${ABACUSSUMMIT_PERSIST}/AbacusSummit/${MINISUITE}/${SIM_NAME}/abacus.par
#echo "Starting ICs using parameter file: $PARFN"
#echo "Using zeldovich module in: $ABACUS"
#$ABACUS/Abacus/zeldovich.py --no-parallel $PARFN
#echo -e "\n\nIC creation complete."
######################
######################

# Merge into one path
SIM_DIR="$SET_NAME/$SIM_NAME"

# Here we assume there are no other sims running on this node!
RAMDISK="/dev/shm/$USER/"

export ABACUS_TMP=$RAMDISK
export ABACUS_SSD=$RAMDISK
ABACUS_SOURCE=$ABACUS


echo -e "ABACUS_PERSIST set to ${ABACUS_PERSIST}"  
export ABACUS_DERIVS=${PROJWORK} # extra env. var. to search for Derivatives on GPFS and copy to BB.
echo -e "ABACUS_DERIVS set to ${ABACUS_DERIVS}"  

echo $PYTHONPATH
echo $PATH

######################################################################## 
######################################################################## 
########################################################################
export PAR_DIR=${ABACUSSUMMIT_PERSIST}/AbacusSummit/${MINISUITE}/${SIM_NAME} #this is where the lsf script that invokes this script copied the pre-processed params. /ccs/ is read-only. 

echo -e "* Running abacus: ${PAR_DIR}\n"  
$ABACUS/Production/run_sim.py "${PAR_DIR}" abacus.par --clean --maxsteps 4

######################################################################## 
######################################################################## 
######################################################################## 

ABACUS_EXIT=$?
echo -e "Abacus exit code: $ABACUS_EXIT \n"  
echo -e "\n\n\n\n"  

#######################################################################
DISBATCH_PY=$ABACUS/external/disBatch/disBatch.py

LOGDIR=$PROJWORK/$USER/logs/$SIM_NAME
mkdir -p $LOGDIR

export DISBATCH_SSH_NODELIST=$(hostname):4

##### Queue up the time slice processing #####
echo -e "About to queue TS processing:"

SIMDIR=$ABACUS_PERSIST/$SET_NAME/$SIM_NAME
SLICEDIRS=$(cd $SIMDIR; find -maxdepth 1 -name 'slice*.*')

if [[ -n "$SLICEDIRS" ]]; then

    # Write a disBatch file with one line per chunk per file type per time slice
    # Chunks should be specified with a literal pattern like slice1.000/*slab001?.L0_pack9.dat
    #DISBATCH_TASKFILE=${SIMDIR}/tmp/${SIM_NAME}.timeslice.disbatch
    
    mkdir $LOGDIR/tmp/
    DISBATCH_TASKFILE=$LOGDIR/tmp/${SIM_NAME}.timeslice.disbatch

    TSSCRIPT=$ABACUS/Production/AbacusSummit/rhea_post_process_timeslice.py
    echo "#DISBATCH PREFIX cd ${SIMDIR}; ${TSSCRIPT} \"" > $DISBATCH_TASKFILE  # write 

    echo "#DISBATCH SUFFIX \" > ${SIMDIR}/\${DISBATCH_NAMETASKS}_\${DISBATCH_JOBID}_\${DISBATCH_TASKID}.log 2>&1" >> $DISBATCH_TASKFILE  # append
    
   #echo "#DISBATCH SUFFIX \" > $LOGDIR/tmp/test.log 2>&1" >> $DISBATCH_TASKFILE  # append


    for SLICE in $SLICEDIRS; do
        # Write one line for every chunk of 10 for every file type
        # slice1.000/*slab001?.L0_pack9.dat

        i=0
        while : ; do
            FNPAT="$(printf '*slab%03d?*.dat' $i)"
            MATCHES=$(cd $SIMDIR; find $SLICE -name "$FNPAT" -print -quit)
            if [[ -z "$MATCHES" ]]; then
                break
            fi
            for FTYPE in L0_pack9 L0_pack9_pids field_pack9 field_pack9_pids; do
                echo "${SLICE}"/$(printf "*slab%03d?.${FTYPE}.dat" $i) >> $DISBATCH_TASKFILE
            done
            i=$((i + 1))
        done
    done

    mkdir -p $LOGDIR/disbatchTS

    TSEPI_SCRIPT=$ABACUS/Production/AbacusSummit/rhea_post_process_timeslice_epilogue.sh
    
    $DISBATCH_PY -e -p $LOGDIR/disbatchTS/$SIM_NAME $DISBATCH_TASKFILE
    $TSEPI_SCRIPT $SIMDIR/slices && echo "Time slices completed successfully."    
fi
####

echo "Done with timeslices." 


##### queue up group processing #####

# mkdir -p tmp
# DISBATCH_TASKFILE=$(pwd)/tmp/${SIM_NAME}.group.disbatch
# DISBATCH_SSH_NODELIST=$(hostname):4

# CHUNK=50
# NWORKER=5
# GROUPSCRIPT="$ABACUS/Abacus/convert_raw_groups_to_asdf.py --chunk=$CHUNK --nworkers=$NWORKER --delete"
# GROUPDIR=$ABACUSSUMMIT_PERSIST/$SET_NAME/$SIM_NAME/group


# if [ -d $GROUPDIR ]; then
#     echo "#DISBATCH PREFIX cd $GROUPDIR; $GROUPSCRIPT " > $DISBATCH_TASKFILE  # write
#     echo "#DISBATCH SUFFIX > $(pwd)/\${DISBATCH_NAMETASKS}_\${DISBATCH_JOBID}_\${DISBATCH_TASKID}.log 2>&1" >> $DISBATCH_TASKFILE  # append

#     GROUPLIST=$(cd $GROUPDIR; echo Step*/ | tr " " "\n")
#     NGROUPLIST=$(wc -l <<< "$GROUPLIST")
#     echo "$GROUPLIST" >> $DISBATCH_TASKFILE

#     mkdir -p ./logs/$SIM_NAME/disbatchGroups
    
#     $DISBATCH_PY -e -p ./logs/$SIM_NAME/disbatchGroups/$SIM_NAME $DISBATCH_TASKFILE && echo Groups completed successfully.
# fi

# echo -e "Completed groups disbatch." 
#######################################################################

export GPFS_PERSIST=$PROJWORK/$USER/nvme/
echo -e "Copying from NVME to ${GPFS_PERSIST}"  

cp -R ${BBPATH}/${SET_NAME}/${SIM_NAME} ${GPFS_PERSIST}/${SIM_NAME}

echo -e "Run and copy to GPFS complete for box ${SIM_NAME} on rank ${JSM_NAMESPACE_RANK}!"  
