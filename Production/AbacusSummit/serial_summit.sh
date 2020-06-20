#!/bin/bash


set -e

export JOB_ACTION_WARNING_TIME=1380

echo "Logs output to :" $LSB_HOSTS

# This script runs on the summit compute node
# Fail if any command returns a non-zero exit status

#Specify simulation to run on this node.
######################
export MINISUITE="SmallBox"
PHASE=$((${JSM_NAMESPACE_RANK} + 3000))
export SIM_NAME="AbacusSummit_small_c000_ph${PHASE}"

SET_NAME="AbacusSummit"
if [[ -z "$SIM_NAME" ]]; then
    echo "Need $SIM_NAME!" ; exit 1
fi


echo -e "* Loading modules\n"  
source modules.sh
echo -e "\n\n\n\n"  

export BBPATH=/mnt/bb/$USER
export ABACUS_PERSIST=$BBPATH
export ABACUS_DERIVS=${PROJWORK} # extra env. var. to search for Derivatives on GPFS and copy to BB.
# Here we assume there are no other sims running on this node!
RAMDISK="/dev/shm/$USER/"
export ABACUS_TMP=$RAMDISK
export ABACUS_SSD=$RAMDISK
ABACUS_SOURCE=$ABACUS

echo -e "If Abacus dir is different from source, making local copy."
export ABACUS=$BBPATH/abacus
# If the Abacus directory is different from the source, make a local copy
if [[ ! $ABACUS -ef $ABACUS_SOURCE ]]; then 
    echo -e "* Copying files locally and running make in $(dirname $ABACUS):\n"
    mkdir -p $(dirname $ABACUS)
    rsync -a --exclude /Production/AbacusSummit/logs $ABACUS_SOURCE/ $ABACUS 
    #cp -R  $ABACUS_SOURCE/ $ABACUS
    #make -C $ABACUS distclean
    #cd $ABACUS
    #./configure CXX=$CXX
    #cd -
    #XX=$CXX make -C $ABACUS clean all
    echo -e "\n\n\n\n"
fi

LOGDIR=$ABACUS/Production/AbacusSummit/logs/$SIM_NAME
mkdir -p $LOGDIR
DEBUG_FN=$LOGDIR/${SIM_NAME}_production.debug
LOG_FN=$LOGDIR/${SIM_NAME}_production.out
echo "Redirecting output to :" ${LOG_FN}
exec >${LOG_FN} 2>&1


echo -e "* Clearing ramdisk\n"
rm -rf $RAMDISK /dev/shm/lgarrison/ /dev/shm/nmaksimova/ /dev/shm/ast145* > $DEBUG_FN
echo -e "\n\n\n\n"

# These actions will be taken upon exit
function cleanup()
{
    export GPFS_PERSIST=$PROJWORK/$USER/nvme/
    echo -e "CLEANUP: Copying from NVME to ${GPFS_PERSIST}"  
    cp -R ${BBPATH}/${SET_NAME}/${SIM_NAME}/ ${GPFS_PERSIST}
    cp -R ${LOGDIR}/ $GPFS_PERSIST/logs/
    echo -e "CLEANUP: Cleanup copy to GPFS complete for box ${SIM_NAME} on rank ${JSM_NAMESPACE_RANK}!"  

    echo -e "* Environment:\n" >> $DEBUG_FN
    printenv &> $DEBUG_FN
    echo -e "\n\n\n\n" >> $DEBUG_FN

    echo -e "* dmesg:\n" >> $DEBUG_FN
    dmesg -T |tail -n15 >> $DEBUG_FN
    echo -e "\n\n\n\n" >> $DEBUG_FN

    echo -e "* /dev/shm usage:\n" >> $DEBUG_FN
    df -h /dev/shm >> $DEBUG_FN
    echo -e "\n\n\n\n" >> $DEBUG_FN

    echo -e "* Clearing ramdisk\n" >> $DEBUG_FN
    rm -rf $RAMDISK >> $DEBUG_FN
    echo -e "\n\n\n\n" >> $DEBUG_FN
}
trap cleanup EXIT


######################################################################## 
######################################################################## 
########################################################################
export PAR_DIR=${ABACUS}/external/AbacusSummit/Simulations/${MINISUITE}/${SIM_NAME} 
echo -e "* Running abacus: ${PAR_DIR}\n"

$ABACUS/Production/run_sim.py "${PAR_DIR}" --clean

ABACUS_EXIT=$?
echo -e "Abacus exit code: $ABACUS_EXIT \n"  
echo -e "\n\n\n\n"  
#######################################################################
#######################################################################

if [[ !($ABACUS_EXIT -eq 0 || $ABACUS_EXIT -eq 200) ]] ; then
     echo -e "Abacus exit code indicates crash. No post-processing. Saving off fns and exiting."   
     exit 
fi

now=$(date +"%T")
echo "Current time, after abacus exit : $now"


DISBATCH_PY=$ABACUS/external/disBatch/disBatch.py

##### queue up group processing #####

echo "Starting group processing." 

mkdir -p $LOGDIR/tmp/
DISBATCH_TASKFILE=$LOGDIR/tmp/${SIM_NAME}.group.disbatch

export DISBATCH_SSH_NODELIST=$(hostname):4

echo ${DISBATCH_SSH_NODELIST}

CHUNK=405 #only one super slab.
NWORKER=1
GROUPSCRIPT="$ABACUS/Abacus/convert_raw_groups_to_asdf.py --chunk=$CHUNK --nworkers=$NWORKER --delete"
GROUPDIR=$ABACUS_PERSIST/$SET_NAME/$SIM_NAME/group

if [ -d $GROUPDIR ]; then
     echo "#DISBATCH PREFIX cd $GROUPDIR; $GROUPSCRIPT " > $DISBATCH_TASKFILE  # write
     echo "#DISBATCH SUFFIX > $LOGDIR/tmp/${SIM_NAME}.log 2>&1" >> $DISBATCH_TASKFILE  # append
     ####echo "#DISBATCH SUFFIX > $LOGDIR/tmp/\${DISBATCH_NAMETASKS}_\${DISBATCH_JOBID}_\${DISBATCH_TASKID}.log 2>&1" >> $DISBATCH_TASKFILE  # append #NAM this consistently fails in disBatch b/c log fn not found

     GROUPLIST=$(cd $GROUPDIR; echo Step*/ | tr " " "\n")
     NGROUPLIST=$(wc -l <<< "$GROUPLIST")
     echo "$GROUPLIST" >> $DISBATCH_TASKFILE

     mkdir -p $LOGDIR/disbatchGroups
    
     $DISBATCH_PY -e -p $LOGDIR/disbatchGroups/ $DISBATCH_TASKFILE && GROUPS_COMPLETED="True"
fi

now=$(date +"%T")
echo "Current time, after groups pp : $now"

if [ "$GROUPS_COMPLETED" = True ] ; then
	echo -e "Completed groups disbatch." 

	WORKINGDIR=$ABACUS_PERSIST/$SET_NAME/$SIM_NAME
	pushd $WORKINGDIR > /dev/null  # work inside sim dir

	echo "Executing epilogue on $WORKINGDIR"

	# Zip and delete the logs in groups of 100 steps
	echo "Starting zip of logs at $(date)"
	pushd log/ > /dev/null
	shopt -s nullglob

	# Delete symlink
	rm -f last

	i=0
	while : ; do
	    start=$(printf 'step%02d??' $i)
	    LOGSTEPDIRS=($start/)
	    if [ "${#LOGSTEPDIRS[@]}" -eq 0 ]; then
	        break
	    fi
	    name=$(printf 'log%02d00' $i)
	    zip --move -qr ${name}.zip ${LOGSTEPDIRS[@]} &
	    i=$((i + 1))
	done

	wait  # wait for all zip to finish
	shopt -u nullglob
	echo "Done zip at $(date)"

	FAST_CKSUM=$ABACUS/external/fast-cksum/bin/fast_cksum

	# Checksum the logs
	$FAST_CKSUM * > checksums.crc32


	now=$(date +"%T")
echo "Current time, after log zip : $now"

	popd > /dev/null  # back to WORKINGDIR

	# Now delete:
	# - Wisdom file: no sense in carting around
	# - Read directory (?): the final header is in the halos anyway, and the full state is in the logs
	# - Any disbatch files, or maybe we'll want to be able to check on them later
	# - Retrieved state
	# - Core files
	echo "Beginning deletions at $(date)"

	echo "Deleting read"
	rm -rf read/

	echo "Deleting core files"
	rm -f core.*

	echo "Deleting retrieved state"
	rm -rf ../${SIM_NAME}_retrieved_state

	echo "Checksum info"
	pushd info/ > /dev/null
	$FAST_CKSUM * > checksums.crc32
	popd > /dev/null

	echo "Deleting extra files"
	rm -f fftw_*.wisdom

	$FAST_CKSUM abacus.par status.log > checksums.crc32

	popd > /dev/null  # back to submission dir
	rm -f tmp/${SIM_NAME}.*
	rm -f slurm_job_*_resize.{sh,csh}

	echo "Done deletions at $(date)"

	# This is a special string the assembly line looks for to know that it's safe to start Globus and htar
	echo "Post processing complete."
        echo -e "SUCCESS: Run complete for box ${SIM_NAME} on rank ${JSM_NAMESPACE_RANK}!"

now=$(date +"%T")
echo "Current time,end : $now"

else
	echo "Groups post processing failed. Cleaning up and exiting."

fi
