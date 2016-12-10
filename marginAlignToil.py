#!/usr/bin/env python
from __future__ import print_function
import os
import uuid
from argparse import ArgumentParser
from toil.job import Job
from toil.common import Toil
import toil_lib.programs as tlp


def cullFilesToLocal(job, fileIds_to_get):
    # type: (toil.job.Job, list<string>) 
    """Gets all the files in 'fileIds_to_get' from the global fileStore
    and puts them in the local working directory. 
    returns the path to the working directory to be passed to docker -v 
    """

    file_dict = {}
    for f in fileIds_to_get:
        file_dict[f] = ""
    
    #work_dir = job.fileStore.getLocalTempDir()
    work_dir = os.path.dirname(os.path.realpath(__file__))
    job.fileStore.logToMaster("[cullDockerFiles] Got directory {}".format(work_dir))
    for file_id in file_dict.keys():
        job.fileStore.logToMaster("preparing {}".format(file_id))
        uniq_file_name = uuid.uuid4().hex + ".tmp"
        destination_path = work_dir + "/{fn}".format(fn=uniq_file_name)
        job.fileStore.readGlobalFile(file_id, userPath=destination_path)
        assert(os.path.exists(destination_path))
        file_dict[file_id] = (destination_path, uniq_file_name)

    return file_dict, work_dir


def bwa_docker_call(job, reference_file_id, reads_file_id, debug=True, memory="10M", cores=1, disk="10M"):
    # type: (toil.job.Job, string, string)
    # get a local copy of the reference file for docker
    file_path_dict, work_dir = cullFilesToLocal(job, [reference_file_id, reads_file_id])
    assert(len(file_path_dict) == 2)
    
    if debug:
        for k in file_path_dict:
            job.fileStore.logToMaster("Looking for temp of {id} at {loc}".format(id=k, loc=file_path_dict[k][0]))
            assert(os.path.exists(file_path_dict[k][0]))
    
    # get the arguments for the docker run command 
    # reference and read file names
    docker_dir = "/data/"
    dkr_reference_path   = docker_dir + file_path_dict[reference_file_id][1]
    bwa_index_parameters = "index " + dkr_reference_path

    bwa_docker_image = "quay.io/ucsc_cgl/bwa"
    
    job.fileStore.logToMaster("workDir: {}".format(work_dir))
    tlp.docker_call(tool=bwa_docker_image, 
                    parameters=bwa_index_parameters.split(), 
                    work_dir=work_dir)

    return


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
            ref_import_string = "file://{abs_path}".format(abs_path=os.path.abspath(args.reference))
            query_import_string = "file://{abs_path}".format(abs_path=os.path.abspath(args.reads))
            print("{r}\n{q}".format(r=ref_import_string, q=query_import_string))
            reference_file_id = toil.importFile(ref_import_string)
            query_file_id = toil.importFile(query_import_string)
            print("{r}\n{q}".format(r=reference_file_id, q=query_file_id))
            root_job = Job.wrapJobFn(bwa_docker_call, reference_file_id, query_file_id)
            return toil.start(root_job)
        else:
            toil.restart()

if __name__ == "__main__":
    main()
