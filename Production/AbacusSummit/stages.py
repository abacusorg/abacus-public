'''
Defines the stages in the assembly line.
The registry of stages is defined in assembly_line.py;
this keeps the flow local to that file and allows us
to define it at the top of the file as well.

While we could trivially move the text descriptions
of the stages from assembly_line.py to each class
below, retaining the text as part of that dictionary
makes it more readable

In general, we probably want to try to catch errors
rather than raise them.  If we raise errors anywhere,
the entire script comes to a halt and we don't write
out the CSV with updated stages or send any emails.
In some cases, it may be appropriate to place the
box in an error stage so the user has to check on it;
in others, it may be appropriate to just not advance
the stage so the script can try again later.
'''

__all__ = ['Stage', 'ErrorStage', 'Frozen',
            'NotStarted', 'ReadyForIC', 'QueuedIC',
            'ReadyForSummit', 'QueuedOnSummit', 'ReadyForPostProcess',
            'QueuedForPostProcess', 'ReadyForDataTransfer', 'QueuedForDataTransfer',
            'ReadyForDeletion','Completed',
            'QueuedJobs', 'globus_status_log']

import os
import os.path
from os.path import join as pjoin, dirname, abspath
import shlex
import subprocess
import socket
import time
import re
import shutil

from Abacus import globus

globus_status_log = 'globus_status.toml'
globus_status_log = pjoin(dirname(abspath(__file__)), globus_status_log)

# Common directory to dump all job script output
logdir = pjoin(dirname(__file__), 'logs')

class Stage:
    # stages are supposed to be immutable, but let's not go overboard enforcing that
    #__slots__ = []

    def __init__(self, name, disable_automation=False, noop=False):
        self.name = name
        self.disable_automation = disable_automation
        self.noop = noop

        # By default, stages have no action, i.e. are passive indicator stages
        if not hasattr(self, 'action'):
            self.action = None

        # One can no-op the action to explicitly disable it
        # noop can be 'action' to no-op the action, or True to no-op the indicator too
        if noop is True:
            self.action = lambda s,b: False  # returns False; stages will not advance
            self.indicator = lambda b: False  # later stages that are no-op'd should not grab boxes
        elif noop == 'action':
            self.action = None

    # Default indicator is no-op
    # Indicators return True if this stage has already been executed, else False
    def indicator(self, box):
        return False

class ErrorStage(Stage):
    def action(self, box):
        '''Don't move on to the next stage if in ErrorStage'''
        return False

class Frozen(Stage):
    def action(self, box):
        '''Don't move on to the next stage if Frozen'''
        return False

#########################################################################################################
################################### STAGE PRECONDITIONS AND ACTIONS #####################################

################################### NOT STARTED ###################################

class NotStarted(Stage):
    def indicator(self, box):
        '''Is OutputDirectory missing, or does it only contain the parameter files?'''

        if self.disable_automation:
            breakpoint()
        try:
            return set(os.listdir(box.params['OutputDirectory'])) == {'abacus.par2', 'abacus.par'}
        except FileNotFoundError:
            return True

    def action(self, box):
        '''reserve IC space'''

        # num_ics_queued   = len([i for i in QueuedJobs('rhea') if 'ICs' in i])

        # size_ics_on_disk = 0
        # for dirpath, dirnames, filenames in os.walk(box.suite_prod_dir):
        #     for f in filenames:
        #         fp = pjoin(dirpath, f)
        #         # tally up ic file sizes; skip if this file is a symbolic link
        #         if not os.path.islink(fp) and f.startswith('ic'):
        #             size_ics_on_disk += os.path.getsize(fp)
        # ic_space_in_use = size_ics_on_disk + IC_DISK_PER_BOX * num_ics_queued
        # # Subtract the total from the budgeted IC disk space. Is there enough_ic_space for this box?
        # if (IC_DISK_RESERVED - ic_space_in_use) > IC_DISK_PER_BOX:

        # TODO: where to get prod_dir?
        num_boxes_ic = 0
        #for direc in [f.path for f in os.scandir(box.suite.suite_prod_dir) if f.is_dir()]:
        #    num_boxes_ic += 'ic' in [f.name for f in os.scandir(direc) if f.is_dir()]

        if num_boxes_ic < 8:
            if self.disable_automation:
                breakpoint()
            os.makedirs(box.params['InitialConditionsDirectory'], exist_ok=False)
            return True
        return False


################################### READY FOR IC ###################################

class ReadyForIC(Stage):
    def indicator(self, box):
        '''Does the IC directory exist?'''
        if self.disable_automation:
            breakpoint()
        return os.path.exists(box.params['InitialConditionsDirectory'])

    def action(self, box):
        '''queue IC creation script on rhea'''
        os.makedirs(f'logs/{box.name}', exist_ok=True)
        cmd = 'sbatch --job-name=' + box.jobname('ICs') + f' -o logs/{box.name}/%x.out rhea_ic_creation.slurm ' + box.parfn
        if self.disable_automation:
            breakpoint()
        try:
            subprocess.run(shlex.split(cmd), check=True)
        except:
            print(f'Error submitting box {box.name} for ICs.  Continuing...')
            return False
        return True

################################### QUEUED FOR IC ###################################

class QueuedIC(Stage):
    def indicator(self, box):
        '''is the IC creation script queued on rhea?'''
        if self.disable_automation:
            breakpoint()
        return box.jobname('ICs') in QueuedJobs('rhea')

################################### READY FOR SUMMIT ###################################

class ReadyForSummit(Stage):
    def indicator(self, box):
        '''are the ICs finished?'''
        if self.disable_automation:
            breakpoint()
        try:
            with open(pjoin(logdir, box.name, box.jobname('ICs') + '.out')) as f:
                return 'IC creation complete.' in f.read()
                #TODO could check ic files sizes, are they what we expect?
        except FileNotFoundError:
            return False

    def action(self, box):
        '''send a signal to summit sleeper script to queue abacus on summit!'''
        success = SignalSummitSleeper('submit', box)
        if not success:
            # box.place_in_error_stage()  # do we need to error out here? maybe not
            return False
        return True

################################### QUEUED ON SUMMIT ###################################

class QueuedOnSummit(Stage):
    def indicator(self, box):
        '''is abacus for this box in the summit queue?'''
        return box.summit_queue_status

        # how does the summit script clean up after itbox to enure it doesn't queue the same box twice?


################################### READY FOR POST-PROCESSING ###################################

class ReadyForPostProcess(Stage):
    def indicator(self, box):
        '''check status.log for "Final Redshift reached, terminating normally"'''
        if self.disable_automation:
            breakpoint()
        try:
            with open(pjoin(box.params['OutputDirectory'], 'status.log')) as f:
                return 'Final redshift' in f.read()
        except FileNotFoundError:
            return False

    def action(self, box):
        '''submit rhea post processing script'''
        cmd = './rhea_post_process.sh ' + box.name
        if self.disable_automation:
            breakpoint()

        try:
            subprocess.run(shlex.split(cmd), check=True)  #TODO slurm python bindings?
        except Exception as e:
            # TODO: do we want to crash here? place box in error stage?
            # Don't really want to raise exception, because then the whole script halts and we don't record anything or send any emails
            print(e)
            print(f'Error submitting box {box.name} for post-procesing.  Continuing...')
            return False

        #also, delete this box's ics! we don't need them anymore.
        return True

################################### QUEUED FOR POST-PROCESSING###################################

class QueuedForPostProcess(Stage):
    def indicator(self, box):
        '''is this box queued for post processing on rhea?'''
        if self.disable_automation:
            breakpoint()
        return box.jobname('PostProcessEpilogue') in QueuedJobs('rhea')

################################### READY FOR DATA TRANSFER ###################################

class ReadyForDataTransfer(Stage):
    def indicator(self, box):
        '''Do we see the "Post processing complete" string in the job output?'''
        if self.disable_automation:
            breakpoint()
        try:
            with open(pjoin(logdir, box.name, box.jobname('PostProcessEpilogue') + '.out')) as f:
                return 'Post processing complete.' in f.read()
        except FileNotFoundError:
            return False

    def action(self, box):
        '''Start Globus and htar transfer'''
        if self.disable_automation:
            breakpoint()

        # the endpoints are actually the defaults, but let's be explicit
        try:
            globus.start_globus_transfer(box.params['OutputDirectory'],
                      dest_path=globus.DEFAULT_NERSC_DEST,
                      source_endpoint=globus.OLCF_DTN_ENDPOINT,
                      dest_endpoint=globus.NERSC_DTN_ENDPOINT,
                      status_log_fn=globus_status_log)
        except Exception as e:
            print(e)
            print(f'Error starting Globus transfer for {box.name}.  Continuing...')
            return False

        # Start an htar job on rhea
        # TODO: for maximum robustness, this should be a separate stage that runs in parallel with the globus stage
        cmd = f'sbatch -M dtn -o logs/{box.name}/' + box.jobname('htar') + '.out --job-name=' + box.jobname('htar') + ' rhea_htar.slurm ' + box.name
        try:
            subprocess.run(shlex.split(cmd), check=True)
        except:
            print(f'Error submitting box {box.name} for htar.  Continuing...')
            return False

        return True

################################### QUEUED FOR DATA TRANSFER ###################################

class QueuedForDataTransfer(Stage):
    def indicator(self, box):
        '''Is Globus or htar queued or running?'''
        gstat = globus.status(box.name, status_log_fn=globus_status_log)
        if gstat is not None and gstat not in globus.GLOBUS_COMPLETION_STATUSES:
            return True

        if box.jobname('htar') in QueuedJobs('dtn'):
            return True

        return False

################################### READY FOR DELETION ###################################

class ReadyForDeletion(Stage):
    def indicator(self, box):
        '''Check that Globus and htar transfers completed successfully'''
        if self.disable_automation:
            breakpoint()

        if globus.status(box.name, status_log_fn=globus_status_log) != 'SUCCEEDED':
            return False

        try:
            with open(pjoin(logdir, box.name, box.jobname('htar') + '.out')) as f:
                if 'Abacus htar complete.' not in f.read():
                    return False
        except FileNotFoundError:
            return False

        return True


    DELETE = ['halos/', 'lightcones/']  # log/
    def action(self, box):
        '''cross your fingers and delete some stuff!'''
        for path in self.DELETE:
            path = pjoin(box.params['OutputDirectory'], path)
            if self.disable_automation:
                breakpoint()
            try:
                shutil.rmtree(path)
            except FileNotFoundError:
                pass  # Not all sims have time slices, e.g.

        # Now touch a file to let us know the lack of files is because we're done this sim
        with open(pjoin(box.params['OutputDirectory'], 'COMPLETED'), 'w') as fp:
            fp.write(time.ctime())

        return True


################################### COMPLETED ###################################

class Completed(Stage):
    def indicator(self, box):
        '''check if output directory contains the COMPLETED file'''
        if self.disable_automation:
            breakpoint()

        return os.path.isfile(pjoin(box.params['OutputDirectory'], 'COMPLETED'))


#########################################################################################################
######################################## EXTRA TIDBITS FOR TIDINESS #####################################

#IMPORTANT! WE RELY ON ALL JOB NAMES BEING UNIQUE. ADD STRINGENT CHECKS FOR THIS.
#TODO: is a queued box ever momentarily not in the queue, e.g. right after submission? or between requeues?
def QueuedJobs(computer_name):
    assert computer_name in ('rhea','summit','dtn')  #sanity check

    # On both summit and rhea, for now we probably only need to know if a job is in queue,
    # don't care about state otherwise.
    # So we'll just return a list of jobs, not a dict

    if computer_name in ('rhea', 'dtn'):
        queued_cmd  = f'squeue -M {computer_name} -A AST145 --format=%j -h'  # -h: no header
        try:
            queued_jobs = subprocess.run(shlex.split(queued_cmd), check=True, capture_output=True, text=True).stdout
        except:
            print(f'Error getting queued jobs on {computer_name}.  Continuing...')
            return []
        return [x.strip() for x in queued_jobs.split('\n')[1:]]  # skip first line, has cluster name

    elif computer_name == 'summit':
        response = SignalSummitSleeper('check')
        if response == ['No unfinished job found']:
            return []

        queued_jobs = [l.split()[0] for l in response]
        return queued_jobs


SUMMIT_SLEEPER_MAGIC_VARS = {'reqfn': 'SUMMIT_SLEEPER_REQUEST',
                             'respfn': 'SUMMIT_SLEEPER_RESPONSE',
                             'submit': 'SubmitToQueue',
                             'check': 'CheckQueue'}
def SignalSummitSleeper(task, box=None, time_limit=300):
    reqfn = SUMMIT_SLEEPER_MAGIC_VARS['reqfn']
    respfn = SUMMIT_SLEEPER_MAGIC_VARS['respfn']
    tmpfn = reqfn + '.tmp'

    taskstr = SUMMIT_SLEEPER_MAGIC_VARS[task]
    if task == 'submit':  # TODO: could aggregate submissions
        taskstr += ':' + box.name

    # Clear any stale responses
    try:
        os.remove(respfn)
    except FileNotFoundError:
        pass

    # Write, flush, sync, close, then rename in hopes of achieving more atomicity
    with open(tmpfn, "w") as f:
        #send a request to the summit sleeper script.
        f.write(taskstr)
        f.flush()
        os.fsync(f.fileno())
    os.rename(tmpfn,reqfn)

    print(f'Sent "{task}" request to Summit sleeper, waiting for response...', flush=True)

    starttime = time.time()

    while time.time() - starttime < time_limit:  #wait for summit sleeper script to return an answer
        try:
            with open(respfn, "r") as f:
                response = f.readlines()
        except FileNotFoundError:
            time.sleep(0.5)
            continue
        os.remove(respfn)

        # sleeper is done!
        if task == 'check':
            # Return the whole string, will be parsed by QueuedJobs
            return response
        elif task == 'submit':
            # Return True or False so the upstream action knows whether to proceed to the next stage
            return bool(re.match(r'Job <\d+> is submitted', response[0]))

    # Timed out waiting for a response
    return False
