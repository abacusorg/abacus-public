#!/usr/bin/env bash

set -e  # Exit on error

# Write stdout and stderr to $LOGFILE
LOGFILE=${0/sh/log}
echo "Sleeper started, logging to $LOGFILE"
exec > $LOGFILE 2>&1
echo "[$(date)] Sleeper started in $(pwd) on host $(hostname)"

REQFN='SUMMIT_SLEEPER_REQUEST'
RESPFN='SUMMIT_SLEEPER_RESPONSE'

submit_job_str='SubmitToQueue'
check_queue_str='CheckQueue'

while :; do  # loop forever
    if [[ -e $REQFN ]]; then
        echo "[$(date)] Waking up..."
        COMMAND=$(<$REQFN)
        rm -f $REQFN

        echo "Received command \"$COMMAND\""
        
        if [[ $COMMAND == $submit_job_str* ]]; then

            # extract the box from the command
            IFS=: read TASK BOX <<< $COMMAND

            echo "Submitting box $BOX to Summit queue."
            # Unlike slurm, lsf doesn't accept command-line script arguments
            RESPONSE=$(SIM_NAME=$BOX bsub -J $BOX -o logs/${BOX}_production.out production.lsf 2>&1)
        
        elif [[ $COMMAND == $check_queue_str* ]]; then
            echo "Checking Summit queue"
            RESPONSE=$(bjobs -o "NAME STAT" -noheader 2>&1)
        else
            # Presumably we want to hard-fail here?
            echo "Unknown task \"$TASK\"!"
            break
        fi

        # Respond to the assembly line
        echo "Responding with: \"$RESPONSE\""
        echo "$RESPONSE" > $RESPFN.tmp
        sync
        mv $RESPFN.tmp $RESPFN

        echo "[$(date)] Done, going back to sleep"
    fi
    # Check again in X seconds
    sleep 5
done

echo "[$(date)] Quitting..."
