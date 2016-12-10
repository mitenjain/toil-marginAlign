#!/usr/bin/env python
from __future__ import print_function
import os
import uuid
import re
import urlparse
from argparse import ArgumentParser
from toil.job import Job
from toil.common import Toil
import toil_lib.programs as tlp


def removeTempSuffix(filename):
    # type: (string)
    """removes '.tmp' from a filename string
    """
    assert filename.endswith(".tmp")
    return re.sub('\.tmp$', '', filename)


def cullFilesToLocal(job, fileIds_to_get):
    # type: (toil.job.Job, list<string>) 
    """Gets all the files in 'fileIds_to_get' from the global fileStore
    and puts them in the local working directory. 
    returns: file_dict<file_id, (full_path, unique_file_name)>, work_dir
    """
    file_dict = {}
    for f in fileIds_to_get:
        file_dict[f] = ""
    
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


def bwa_index_docker_call(job, reference_file_id,
                          debug=True, 
                          memory="10M", cores=1, disk="10M", 
                          bwa_docker_image="quay.io/ucsc_cgl/bwa"):
    # type: (toil.job.Job, string, string, 
    #        bool, 
    #        string, int, string, 
    #        string)
    """Builds all required indices (crufty files) that BWA needs and 
    imports them into the fileStore
    """
    # function variables
    DockerDir = "/data/"
    BwaReferenceFiles = [reference_file_id]
    
    def _get_unique_filename(file_id):
        return file_path_dict[file_id][1]
    
    def _get_unique_filepath(file_id):
        return file_path_dict[file_id][0]

    def _safe_delete(local_file_list):
        removed_bools = []
        paths = [urlparse.urlparse(f).path for f in local_file_list]
        job.fileStore.logToMaster("got {fl} turned into {pts}".format(fl=local_file_list, pts=paths))
        for url in paths:
            assert os.path.exists(url), "Didn'``t find {}".format(url)
            os.remove(url)
            removed_bools.append(not os.path.exists(url))
        return removed_bools
    
    def _run_bwa_index():
        bwa_index_parameters = ["index", dkr_reference_path]
        if debug:
            job.fileStore.logToMaster("workDir: {}".format(work_dir))
        # bwa docker call creates index files in the local working directory
        tlp.docker_call(tool=bwa_docker_image, 
                        parameters=bwa_index_parameters,
                        work_dir=work_dir)
        # import files to fileStore
        bwa_index_urls = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        bwa_index_urls = ["file://" + _get_unique_filepath(reference_file_id) + x for x in bwa_index_urls]
        new_ids        = [job.fileStore.importFile(x) for x in bwa_index_urls]
        removed        = _safe_delete(bwa_index_urls)

        if debug:
            for i, x in enumerate(removed):
                if x is False:
                    job.fileStore.logToMaster("[bwa_docker_call::_run_bwa_index] " 
                                              "failed to remove {}".format(bwa_index_urls[i]))
                else:
                    job.fileStore.logToMaster("[bwa_docker_call::_run_bwa_index] "
                                              "removed {}".format(bwa_index_urls[i]))
                    
        return new_ids

    # get a local copy of the reference and reads files for docker
    file_path_dict, work_dir = cullFilesToLocal(job, BwaReferenceFiles)
    assert(len(file_path_dict) == 1)
    
    if debug:
        for k in file_path_dict:
            job.fileStore.logToMaster("Looking for temp of {id} at {loc}".format(id=k, loc=file_path_dict[k][0]))
            assert(os.path.exists(file_path_dict[k][0]))
    
    # arguments for the bwa indexing and alignment
    dkr_reference_path   = DockerDir + _get_unique_filename(reference_file_id)
    BwaReferenceFiles += _run_bwa_index()
    return BwaReferenceFiles


def bwa_docker_align(job, reference_id_pack, reads_file_id, 
                     bwa_docker_image="quay.io/ucsc_cgl/bwa"):
    job.fileStore.logToMaster("BWA reference id pack {}".format(reference_id_pack.__str__())) 



def bwa_docker_alignment_root(job, reference_file_id, reads_file_id, 
                              debug=True, 
                              memory="10M", cores=1, disk="10M", 
                              bwa_docker_image="quay.io/ucsc_cgl/bwa"):
    BwaRefIdPack = job.addChildJobFn(bwa_index_docker_call, reference_file_id).rv()
    job.addFollowOnJobFn(bwa_docker_align, BwaRefIdPack, reads_file_id)
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
            ref_import_string   = "file://{abs_path}".format(abs_path=os.path.abspath(args.reference))
            query_import_string = "file://{abs_path}".format(abs_path=os.path.abspath(args.reads))
            reference_file_id   = toil.importFile(ref_import_string)
            query_file_id       = toil.importFile(query_import_string)
            root_job            = Job.wrapJobFn(bwa_docker_alignment_root, reference_file_id, query_file_id)
            return toil.start(root_job)
        else:
            toil.restart()

if __name__ == "__main__":
    main()
