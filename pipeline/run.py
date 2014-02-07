#!/usr/bin/env python

import sys
import argparse
import pprint

from cosmos import session
from cosmos.Workflow.models import Workflow
import cosmos.Workflow.cli

from orchestra import make_drmaa_lsf_spec
from workflow import build_dag

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    cosmos.Workflow.cli.add_workflow_args(parser)
    wf_args, args = cosmos.Workflow.cli.parse_args(parser)

    dag = build_dag(args['input'], chrs=args['chr'])

    # Run things
    wf = Workflow.start(**wf_args)
    wf.log.info('----- CONFIG -----\n' + pprint.pformat(args))
    session.get_drmaa_native_specification = make_drmaa_lsf_spec(args['name'])
    dag.configure(args)
    dag.add_to_workflow(wf)
    wf.run()
