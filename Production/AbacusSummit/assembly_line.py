#!/usr/bin/env python3

'''
######################## AbacusSummit #######################
#################### assembly line utility ##################

# This script's purpose is to detect the current stages of
# the AbacusSummit boxes and generate a csv file describing
# them. It will then execute next eligible production steps
# automatically. If this script encounters an unusual state,
# it will email Nina & Lehman.

#############################################################
'''

import argparse
import csv
import os
from os.path import join as pjoin, abspath
import smtplib
from email.message import EmailMessage
import shlex
import time

from Abacus import abacus, GenParam, Tools
from Abacus.InputFile import InputFile
from Abacus import globus

# import * often not a good idea, but we defined __all__ in stages.py for this purpose
from stages import *

os.environ['ABACUS_PERSIST'] = os.environ['ABACUSSUMMIT_PERSIST']

#########################################################################################################
############################################## SUITE STAGES #############################################

stages = {
       -2:                Frozen('Frozen by user'),
       -1:            ErrorStage('Error!'),
        0:            NotStarted('Not Started'),
        1:            ReadyForIC('Ready for IC'),
        2:              QueuedIC('Queued IC'),
        3:        ReadyForSummit('Ready for Summit'),
        4:        QueuedOnSummit('Queued on Summit'),
        5:   ReadyForPostProcess('Ready for post-processing'),
        6:  QueuedForPostProcess('Queued for post-processing'),
        7:  ReadyForDataTransfer('Ready for data transfer'),
        8: QueuedForDataTransfer('Queued for data transfer'),
        9:      ReadyForDeletion('Ready for deletion'),
        10:            Completed('Completed'),
    }

# Inform stages of their numbers
for stagenum in stages:
    stages[stagenum].num = stagenum

##TODO: HOW DO WE PREVENT ACTIONS HAPPENING TWICE DESTRUCTIVELY?

class Box:
    def __init__(self, box_dir, suite, disable_automation, summit_queue_status):
        self.par2 = pjoin(box_dir, 'abacus.par2')
        self.disable_automation = disable_automation
        self.parfn = pjoin(box_dir, 'abacus.par')
        self.suite = suite
        self.known_error = False
        self.summit_queue_status = summit_queue_status

        try:  # if abacus.par already exists, import it.
            self.params = InputFile(self.parfn)
        except FileNotFoundError: # otherwise, generate it.
            self.params = abacus.preprocess_params(pjoin(box_dir, 'abacus.par'), self.par2)

        self.name   = self.params['SimName']
        self.stage  = self.detect_stage_num() #this returns a Stage object if valid stage detected, None otherwise.


    def detect_stage_num(self, disable_automation=None):
        if disable_automation is None:
            disable_automation = self.disable_automation

        #first, check if this box has been frozen by the user. This disables all attempts at actions, and is indicated by the prence of a FREEZE file.
        if os.path.exists(pjoin(self.params['OutputDirectory'], 'FREEZE')):
            return stages[-2] #return Frozen stage.

        # by iterating through the possible stages in reverse, the indicators can be a little less stringent
        for stage_num in reversed(sorted(stages)):
            thisstage = stages[stage_num]
            # this relies on indicators being exclusive! but maybe if we iterate in reverse, it's easier.
            if thisstage.indicator(self):
                if self.disable_automation: breakpoint()
                return stages[stage_num]
        return stages[-1]


    def jobname(self, jobtype):
        return self.name + '_' +jobtype

    def attempt_action(self, verbose=True):
        attempted_action = self.stage.num
        if self.stage.action is None:
            if verbose:
                print(f'[{self.name}] Stage "{self.stage.name}" has no action. Continuing.')
                return None

        if verbose:
            print(f'Attempting "{self.stage.name}" action for box {self.name}')
        success = self.stage.action(self)
        if success:
            print("\tAction succeeded.")
            self.stage = stages[self.stage.num + 1]
            return attempted_action
        else:
            print("\tAction did not execute.")
            return None

    # allow comparison of Box objects by box name
    def __eq__(self, other):
        if isinstance(other, str):
            return self.name == other
        elif isinstance(other, Box):
            return self.name == other.name
        return False

    def __ne__(self,other):
        return not (self == other)


class Suite:
    #detect stages of all boxes in suite. generate csv file summarizing them.
    def __init__(self, csv, suite_spec_dir, suite_name='AbacusSummit', add_boxes=None, disable_automation=True,
                    refresh=False, no_summit=False, verbose=True, email=False, no_globus=False):
        self.suite_spec_dir = suite_spec_dir
        self.disable_automation = disable_automation
        self.verbose = verbose
        self.suite_name = suite_name
        self.refresh = refresh
        self.no_summit = no_summit
        self.no_globus = no_globus
        self.email = email

        self.csv_fn = csv
        prev = self.load_previous_status()

        self.summit_queue_status = QueuedJobs('summit') if not no_summit else []

        if not no_globus:
            globus.update_all_status_from_globus(status_log_fn=globus_status_log)

        # By default, only process boxes that are present in the CSV
        self.boxes = []
        if prev is not None:
            for boxname in prev:
                box = Box(pjoin(suite_spec_dir, boxname), self, self.disable_automation, boxname in self.summit_queue_status)
                self.boxes += [box]

        if add_boxes is not None:
            assert len(set(add_boxes)) == len(add_boxes)  # no duplicates in new
            assert all(ab not in [box.name for box in self.boxes] for ab in add_boxes)  # no duplicates with existing
            for boxname in add_boxes:
                box = Box(pjoin(suite_spec_dir, boxname), self, self.disable_automation, boxname in self.summit_queue_status)
                self.boxes += [box]


        #sort boxes by descending stage number, so that we attempt later actions first!
        self.boxes.sort(key=lambda b: b.stage.num, reverse=True)
        self.actions_taken = {}  # {box.name:stage.num}

        self.sanity_check()


    def load_previous_status(self):
        '''Read the CSV file'''
        try:
            with open(self.csv_fn, 'r', newline='') as file:
                reader = csv.reader(file)
                firstline = next(reader)[0]
                # Just to prevent later overwriting the wrong file
                if firstline != '# AbacusSummit CSV File':
                    raise RuntimeError(f"File {self.csv_fn} does not look like an AbacusSummit CSV file")

                #TODO this is dangerous if we ever change record_summary_csv()
                self.previous_suite_status = {row[0]   : int(row[1])   for row in reader if not row[0].startswith('#')}
        except FileNotFoundError:
            self.previous_suite_status = None

        return self.previous_suite_status


    def sanity_check(self):
        if self.previous_suite_status is None:
            return

        current_suite_status  = {box.name : box.stage.num for box in self.boxes }

        for box in self.boxes:
            if self.disable_automation:
                breakpoint()
            this_stage_num = box.stage.num
            # Get the previous stage num from the CSV
            # but if the box is new then it won't be in the CSV so use the current stage
            last_stage_num = self.previous_suite_status.get(box.name, this_stage_num)

            # If the box was placed in an error state last run, don't let it escape that state unless run with --refresh
            # TODO: do we ever want to allow boxes to escape the error state on their own? To overcome transient errors?
            # In that case, we'd need to allow jumps of > 1 when coming out of the error state
            # TODO: do we need to write actions taken to the CSV immediately? If a later box causes something to crash before
            # the CSV is written, then next time we start, we may see unexpected jumps
            if last_stage_num == -1 and not self.refresh:
                box.stage = stages[-1]
                box.known_error = True
                continue

            if box.known_error:
                continue

            subject = f'Box {box.name} needs attention!'

            # TODO: do we want to attempt actions even if some box is in an invalid state?
            if this_stage_num == -1:
                if self.disable_automation: breakpoint()
                self.send_email(subject,  f'No valid stage detected for box: {box.name}\n')
                box.known_error = True
                continue

            if this_stage_num == -2:
                if self.disable_automation: breakpoint()
                continue

            if this_stage_num < last_stage_num:
                if self.refresh and last_stage_num < 0: # only allow boxes in error states to be refreshed. 
                    # TODO: append all(?) print statements to emails and log it
                    print(f'Refresh is allowing box {box.name} to regress from stage {last_stage_num} to {this_stage_num}')
                    continue

                if self.disable_automation: breakpoint()
                self.send_email(subject, f'Error: box {box.name} has regressed from stage {last_stage_num} to {this_stage_num}!')
                box.stage = stages[-1]
                box.known_error = True
                continue

            if this_stage_num > last_stage_num + 2:
                if self.refresh and last_stage_num < 0:
                    print(f'Refresh is allowing box {box.name} to increment from {last_stage_num} to {this_stage_num} by more than two actions')
                    continue
                if self.disable_automation: breakpoint()
                self.send_email(subject, f'Error: box {box.name} has incremented its last stage {last_stage_num} to {this_stage_num} by more than two actions!')
                box.stage = stages[-1]
                box.known_error = True
                continue

            if this_stage_num == last_stage_num + 2: # the only boxes that could have incremented their stage by 2 are the ones that were
                if self.disable_automation: breakpoint()
                last_stage = stages[last_stage_num] # previously at a stage where the action was a NoOp.
                if last_stage.action is not None:
                    if self.refresh and last_stage_num < 0:
                        print(f'Refresh is allowing box {box.name} to increment from {last_stage_num} to {this_stage_num}')
                        continue
                    self.send_email(subject, f'Error: box {box.name} went from stage {last_stage_num} to {this_stage_num}, even though its action is not a NoOp!')
                    box.stage = stages[-1]
                    box.known_error = True
                    continue


    def attempt_actions(self):
        for box in self.boxes:
            if self.disable_automation: breakpoint()
            taken_action = box.attempt_action()
            if taken_action is not None:
                assert box.name not in self.actions_taken
                self.actions_taken[box.name] = taken_action

        self.sanity_check()


    def record_summary_csv(self):
        header = f'# AbacusSummit CSV File\n' + \
                 f'# Suite specification directory: {self.suite_spec_dir}\n' + \
                 f'# Suite working directory: ' + os.environ['ABACUS_PERSIST'] + f'/{self.suite_name}\n'
        subject = '[AbacusSummit] Status'
        message = header
        with open(self.csv_fn, 'w', newline='') as fp:
            fp.write(header)
            writer = csv.writer(fp)
            for box in self.boxes:
                ret = writer.writerow([box.name, box.stage.num])
                message += f'{box.name}: {box.stage.num} ({box.stage.name})\n'

            # If we only updated a subset of the boxes, pass through any untouched boxes
            if self.previous_suite_status:
                for boxname in self.previous_suite_status:
                    if boxname not in self.boxes:
                        ret = writer.writerow([boxname, self.previous_suite_status[boxname]])
                        message += f'{boxname}: {self.previous_suite_status[boxname]} ({stages[self.previous_suite_status[boxname]].name})\n'

        self.send_email(subject, message)


    def send_email(self, subject, message, dryrun=None):
        if dryrun is None:
            dryrun = not self.email

        message += '\nNew actions taken:\n'
        if self.actions_taken:
            for b,anum in self.actions_taken.items():
                message += f'\t{b}: {anum} ({stages[anum].name}: {stages[anum].action.__doc__})\n'
        else:
            message += '\tNone.\n'

        message += f'\n----------\nThis is an automated email sent by the AbacusSummit assembly line script, ' + \
                    f'running at {abspath(__file__)}.  Message generated at {time.asctime()}.'

        print(f'Email\n\nSubject: {subject}\nMessage\n========\n{message}\n========\n')
        if dryrun:
            return

        msg = EmailMessage()
        msg.set_content(message)
        people = 'lgarrison@flatironinstitute.org', 'nina.maksimova@cfa.harvard.edu', 'deisenstein@cfa.harvard.edu'
        msg['Subject'] = subject
        msg['From'] = 'noreply'
        msg['To'] = people
        msg['reply-to'] = people

        # Send the message via SMTP server.
        s = smtplib.SMTP('localhost')
        s.send_message(msg)
        s.quit()

        print(f'Sent email "{subject}"')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter) #TODO add docstring.
    parser.add_argument('csv', help='Name of the CSV file used to track the status of this set of sims')
    parser.add_argument('--manual' , help="Override automation in assembly line script. Ask user permission interactively before attempting next action.",
        action='store_true')
    parser.add_argument('--verbose', help="Enable verbose outputs", action='store_true', default=True)
    parser.add_argument('--spec-dir', help="The suite specification directory containing the directories with .par2 files for each sim, like 'AbacusSummit/Simulations'",
        default=os.environ.get('ABACUSSUMMIT_SPEC'))
    parser.add_argument('--no-action', help="Just update the status of the CSV, don't take any actions", action='store_true')
    parser.add_argument('--no-summit', help="Assume no jobs are queued on summit", action='store_true')
    parser.add_argument('--no-globus', help="Assume no jobs are running on Globus", action='store_true')
    parser.add_argument('--email', help="Send status and error emails", action='store_true')
    parser.add_argument('--add-boxes', help="Add these box names to the current CSV (or make a new CSV with these boxes)", nargs='+')
    parser.add_argument('--refresh', help="Allow boxes to escape the Error stage", action='store_true')

    args = parser.parse_args()
    args = vars(args)
    no_action = args.pop('no_action')
    args['suite_spec_dir'] = args.pop('spec_dir')
    if not args['suite_spec_dir']:
        raise ValueError('Must specify --spec-dir or $ABACUSSUMMIT_SPEC')
    args['disable_automation'] = args.pop('manual')

    #detect stages and do sanity check.
    suite = Suite(**args, suite_name='AbacusSummit')  #email=not no_action

    #attempt actions for all boxes, then do another sanity check.
    if not no_action:
        suite.attempt_actions()

    #record suite status and email users.
    suite.record_summary_csv()


    #TODO if a box is running, add a file that causes crash if any other job attempts to touch it.
