#!/usr/bin/env python3

'''
This script facilitates Globus transfers of Abacus simulation data.
In particular, for AbacusSummit we want to transfer particle data and
halo catalogs from OLCF to NERSC after each simulation is done.
Both sites already have Globus endpoints, so the job of this script
is to initiate a transfer between the two, saving the transfer ID
number so we can query the transfer status later.

The intention of this script is primarily to do disk-to-disk
transfers rather than invoking HPSS, which is probably best
done separately via htar.

Important: this script assumes that a simulation's name is a unique
identifier.

Usage
=====
$ ./globus.py --help

Based on Globus automation-examples/globus_folder_sync.py.
'''

import sys
import os
import urllib
import argparse

import toml

from globus_sdk import (NativeAppAuthClient, TransferClient,
                        RefreshTokenAuthorizer, TransferData)
from globus_sdk.exc import GlobusAPIError, TransferAPIError

from fair_research_login import NativeClient

from Abacus.Tools import ArgParseFormatter

# source and destination endpoints
OLCF_DTN_ENDPOINT = 'ef1a9560-7ca1-11e5-992c-22000b96db58'
NERSC_DTN_ENDPOINT = '9d6d994a-6d04-11e5-ba46-22000b92c6ec'

# Destination Path -- The directory will be created if it doesn't exist
DEFAULT_NERSC_DEST = '/global/cfs/cdirs/desi/cosmosim/Abacus'

# You will need to register a *Native App* at https://developers.globus.org/
# Your app should include the following:
#     - The scopes should match the SCOPES variable below
#     - Your app's clientid should match the CLIENT_ID var below
#     - "Native App" should be checked
# For more information:
# https://docs.globus.org/api/auth/developer-guide/#register-app
CLIENT_ID = 'c8b96d23-7922-40e0-9344-aceaff7ef99d'
REDIRECT_URI = 'https://auth.globus.org/v2/web/auth-code'
SCOPES = ('openid email profile '
          'urn:globus:auth:scope:transfer.api.globus.org:all')

APP_NAME = 'AbacusSummit Data Transfer'

# The status of jobs will not be re-queried if they were in one of these states last time we checked
GLOBUS_COMPLETION_STATUSES = ['SUCCEEDED', 'FAILED']  # other possibilities are ACTIVE, INACTIVE

# Create the destination folder if it does not already exist
CREATE_DESTINATION_FOLDER = True

DEFAULT_STATUS_FILE = 'globus_status.toml'
DEFAULT_STATUS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), DEFAULT_STATUS_FILE)  # make absolute


def load_status_log(status_log_fn):
    try:
        status_log = toml.load(status_log_fn)
    except FileNotFoundError:
        status_log = {}
    return status_log


def write_status_log(status_log, status_log_fn):
    with open(status_log_fn, 'w') as fp:
        toml.dump(status_log, fp)


def update_all_status_from_globus(transfer_client=None, status_log=None,
                                    status_log_fn=DEFAULT_STATUS_FILE,
                                    source_endpoint=OLCF_DTN_ENDPOINT, dest_endpoint=NERSC_DTN_ENDPOINT):
    '''
    For each unfinished box in the status_log dict,
    query globus to get its status_log and update
    the dict.

    If `transfer_client` and `status_log` can be given
    if those resources are already open.  Otherwise
    they will be opened.
    '''

    own_status_log = False
    if status_log is None:
        status_log = load_status_log(status_log_fn)
        own_status_log = True

    if transfer_client is None:
        transfer_client = setup_transfer_client(source_endpoint, dest_endpoint)

    for boxname in status_log:
        boxstatus = status_log[boxname]
        
        oldstatus = boxstatus.get('status',None)

        if oldstatus in GLOBUS_COMPLETION_STATUSES:
            # We already know the fate of this task, no need to check again
            continue

        submission = boxstatus
        task = transfer_client.get_task(submission['task_id'])
        
        task = task.data  # get dict
        del task['event_link']

        if task['status'] != oldstatus:
            print(f'Task status for {boxname} changed from {oldstatus} to {task["status"]}')
        status_log[boxname] = task

    if own_status_log:
        write_status_log(status_log, status_log_fn)


def status(boxname, status_log_fn=DEFAULT_STATUS_FILE):
    '''Get the log status of a single box'''
    status_log = load_status_log(status_log_fn)

    if boxname in status_log:
        return status_log[boxname]['status']
    return None


def get_client_tokens():
    tokens = None
    client = NativeClient(client_id=CLIENT_ID, app_name=APP_NAME)
    try:
        # if we already have tokens, load and use them
        tokens = client.load_tokens(requested_scopes=SCOPES)
    except:
        pass

    if not tokens:
        # if we need to get tokens, start the Native App authentication process
        # need to specify that we want refresh tokens
        # N.B. the no_local_server is the key option not in the Globus automation-examples
        # that lets us accomplish the login on a remote node.
        tokens = client.login(requested_scopes=SCOPES,
                              refresh_tokens=True, no_local_server=True)
        try:
            client.save_tokens(tokens)
        except:
            pass

    return tokens


def setup_transfer_client(source_endpoint, dest_endpoint):
    tokens = get_client_tokens()
    transfer_tokens = tokens['transfer.api.globus.org']

    authorizer = RefreshTokenAuthorizer(
        transfer_tokens['refresh_token'],
        NativeAppAuthClient(client_id=CLIENT_ID),
        access_token=transfer_tokens['access_token'],
        expires_at=transfer_tokens['expires_at_seconds'])

    transfer_client = TransferClient(authorizer=authorizer)

    try:
        transfer_client.endpoint_autoactivate(source_endpoint)
        transfer_client.endpoint_autoactivate(dest_endpoint)
    except GlobusAPIError as ex:
        if ex.http_status == 401:
            # TODO: when does this happen? Does the user need to delete a local file or just login again?
            sys.exit('Refresh token has expired.')
        else:
            raise ex
    return transfer_client


def check_endpoint_path(transfer_client, endpoint, path):
    """Check the endpoint path exists"""
    try:
        transfer_client.operation_ls(endpoint, path=path)
    except TransferAPIError as tapie:
        sys.exit(f'Failed to query endpoint "{endpoint}": {tapie.message}')


def create_destination_directory(transfer_client, dest_ep, dest_path):
    """Create the destination path if it does not exist"""
    try:
        transfer_client.operation_ls(dest_ep, path=dest_path)
    except TransferAPIError:
        try:
            transfer_client.operation_mkdir(dest_ep, dest_path)
            print('Created directory on destination endpoint: {}'.format(dest_path))
        except TransferAPIError as tapie:
            sys.exit(f'Failed to create directory on destination endpoint: {tapie.message}')


def _submit_transfer(tdata, transfer_client, status_log, status_log_fn, boxname):
    '''
    Submit the tdata TransferData object to the transfer_client

    TODO: would probably make sense to wrap this whole interface in an object
    so we don't have to pass so many params
    '''
    
    assert boxname not in status_log  # should never fail, but double check
    task = transfer_client.submit_transfer(tdata)
    
    task = task.data  # turn into a dict
    del task['task_link']  # don't need this in the log file
    task['status'] = task.pop('code')  # rename for consistency with the task status
    status_log[boxname] = task
    
    # Immediately write the status log, this is our only (local) record
    # Of course, one can always check the Globus web interface too
    write_status_log(status_log, status_log_fn)

    assert task['status'] == 'Accepted'


def start_globus_transfer(source_path, dest_path=DEFAULT_NERSC_DEST,
                          source_endpoint=OLCF_DTN_ENDPOINT, dest_endpoint=NERSC_DTN_ENDPOINT,
                          globus_verify_checksum=True, status_log_fn=DEFAULT_STATUS_FILE):
    '''
    Start a Globus transfer from the `source_path` on the `source_endpoint` to `dest_path`
    on the `dest_endpoint`.

    Note that `source_path` and `dest_path` are the paths as seen from the Globus endpoint,
    which might(??) not be the same as the regular filesystem path.  One can double-check
    this with the Globus web interface.
    '''

    # Determine the box name from the file path
    boxname = os.path.basename(os.path.abspath(source_path))
    assert(os.path.isdir(source_path))

    # Globus transfers a directory's contents, not the directory itself
    # so need to make the containing directory on the destination
    dest_path = os.path.join(dest_path, boxname)

    # Check that this box is not already in the status file
    status_log = load_status_log(status_log_fn)
    if boxname in status_log:
        # Could instead check if status is FAILED (or whatever) and allow reruns in that case
        raise RuntimeError(f'Box {boxname} was already found in {status_log_fn}. Is this a duplicate transfer?')

    # Authenticate with Globus
    transfer_client = setup_transfer_client(source_endpoint, dest_endpoint)

    # Set up the destination directory
    check_endpoint_path(transfer_client, source_endpoint, source_path)
    if CREATE_DESTINATION_FOLDER:
        create_destination_directory(transfer_client, dest_endpoint,
                                     dest_path)
    else:
        check_endpoint_path(transfer_client, dest_endpoint, dest_path)

    tdata = TransferData(
        transfer_client,
        source_endpoint,
        dest_endpoint,
        label=f'{boxname} transfer',
        sync_level=None,  # always overwrite
        verify_checksum=globus_verify_checksum,
        encrypt_data=False,
        preserve_timestamp=True,
        recursive_symlinks='ignore'  # not supported
    )
    # TODO: copy dir instead of contents?
    tdata.add_item(source_path, dest_path, recursive=True)

    # also writes status log
    _submit_transfer(tdata, transfer_client, status_log, status_log_fn, boxname)

    print('Transfer has been started from\n  {}:{}\nto\n  {}:{}'.format(
        source_endpoint,
        source_path,
        dest_endpoint,
        dest_path
    ))
    url_string = 'https://app.globus.org/app/transfer?' + \
        urllib.parse.urlencode({
            'origin_id': source_endpoint,
            'origin_path': source_path,
            'destination_id': dest_endpoint,
            'destination_path': dest_path
        })
    print('Visit the link below to see the changes:\n{}'.format(url_string))

    # While we're here, go ahead and update the status of any pending transfers for other boxes
    update_all_status_from_globus(transfer_client, status_log)
    write_status_log(status_log, status_log_fn)


def subcommand_transfer(**args):
    '''Process the command-line args for the "transfer" subcommand'''

    args['globus_verify_checksum'] = not args.pop('no_globus_checksum')

    for sim_dir in args.pop('simulation-dir'):
        start_globus_transfer(sim_dir, **args)


def subcommand_status(simulation, status_log_fn=DEFAULT_STATUS_FILE, **kwargs):
    # In case we were passed a path, get absolute basename
    boxname = os.path.basename(os.path.abspath(simulation))

    gstat = status(boxname, status_log_fn)
    print('Not found in log' if gstat is None else gstat)


if __name__ == '__main__':
    # TODO: want to show main help alongside subcommand help, or show all help in monolithic block
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=ArgParseFormatter)
    subparsers = parser.add_subparsers(help='Sub-commands, git-style', required=True, dest='subcommand')

    parser.add_argument('--status-log', help='Globus status log filename', default=DEFAULT_STATUS_FILE)
    parser.add_argument('--source-endpoint', help='The Globus source endpoint', default='OLCF_DTN_ENDPOINT')
    parser.add_argument('--dest-endpoint', help='The Globus destination endpoint', default='NERSC_DTN_ENDPOINT')

    tparser = subparsers.add_parser('transfer', help='Transfer one or more boxes via Globus',
        description='Transfer one or more boxes via Globus')
    tparser.add_argument('simulation-dir', help='The simulation OutputDirectory containing the time slices and halo catalogs to transfer.', nargs='+')
    tparser.add_argument('--no-globus-checksum', help='Turn off Globus checksum verification', action='store_true')
    tparser.set_defaults(func=subcommand_transfer)

    uparser = subparsers.add_parser('update', help='Update the status log from Globus',
        description='Update the status log by querying Globus for the status of all pending boxes in the log')
    uparser.set_defaults(func=update_all_status_from_globus)  # just a straight pass-through

    sparser = subparsers.add_parser('status', help='Get the log status of one box',
        description='Read the local status log to get the transfer status of a box (based on last time we checked Globus')
    sparser.add_argument('simulation', help='The simulation dir or box name. Only the basename will be used if a path is given.')
    sparser.set_defaults(func=subcommand_status)
    
    args = parser.parse_args()
    func = args.func
    args = vars(args)

    args['status_log_fn'] = args.pop('status_log')
    # Allow user to pass *_ENDPOINT variable names
    for k in ('source_endpoint','dest_endpoint'):
        if args[k].endswith('_ENDPOINT'):
            args[k] = eval(args[k])

    del args['subcommand'], args['func']

    func(**args)
