#!/usr/bin/env python
from __future__ import print_function
import os
from argparse import ArgumentParser
from toil.job import Job
from toil.common import Toil
from bwa import bwa_docker_alignment_root


def main():
    def parse_args():
        parser = ArgumentParser()
        parser.add_argument("--reference", "-r", dest="reference", required=True)
        parser.add_argument("--reads", "-q", dest="reads", required=True)
        parser.add_argument("--out_sam", "-o", dest="out_sam", required=True)
        Job.Runner.addToilOptions(parser)
        return parser.parse_args()

    args = parse_args()
    
    with Toil(args) as toil:
        if not toil.options.restart:
            ref_import_string    = "file://{abs_path}".format(abs_path=os.path.abspath(args.reference))
            query_import_string  = "file://{abs_path}".format(abs_path=os.path.abspath(args.reads))
            output_export_string = "file://{abs_path}".format(abs_path=os.path.abspath(args.out_sam))
            reference_file_id    = toil.importFile(ref_import_string)
            query_file_id        = toil.importFile(query_import_string)
            root_job             = Job.wrapJobFn(bwa_docker_alignment_root, 
                                                 reference_file_id, query_file_id, 
                                                 output_export_string)
            return toil.start(root_job)
        else:
            toil.restart()

if __name__ == "__main__":
    main()
