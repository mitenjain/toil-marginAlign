#!/usr/bin/env python
from __future__ import print_function
import os
import uuid
import re
import urlparse
from argparse import ArgumentParser
from itertools import izip
from toil.job import Job
from toil.common import Toil
import toil_lib.programs as tlp

DEBUG = True
DOCKER_DIR = "/data/"


def removeTempSuffix(filename):
    # type: (string)
    """removes '.tmp' from a filename string
    """
    assert filename.endswith(".tmp")
    return re.sub('\.tmp$', '', filename)


def bwa_index_file_suffixes():
    return [".amb", ".ann", ".bwt", ".pac", ".sa"]


class LocalFileManager(object):
    """Gets all the files in 'fileIds_to_get' from the global fileStore
    and puts them in the local working directory. 
    returns: file_dict<file_id, (full_path, unique_file_name)>, work_dir
    """
    def __init__(self, job, fileIds_to_get, userFileNames=None):
        # type: (toil.job.Job, list<str>, dict<str, str>)
        self.owner_job = job
        self.work_dir  = os.path.dirname(os.path.realpath(__file__))
        self.file_dict = self._initialize_file_dict(fileIds_to_get, userFileNames)

    def localFileName(self, fileId):
        return self.file_dict[fileId][1]

    def localFilePath(self, fileId):
        return self.file_dict[fileId][0]

    def workDir(self):
        return self.work_dir + "/"

    @staticmethod
    def safeLocalDelete(local_file_list):
        removed_bools = []
        paths = [urlparse.urlparse(f).path for f in local_file_list]
        for url in paths:
            assert os.path.exists(url), "Didn'``t find {}".format(url)
            os.remove(url)
            removed_bools.append(not os.path.exists(url))
        return removed_bools

    def _initialize_file_dict(self, file_ids, userFileNames):
        file_dict = {}
        for f in file_ids:
            file_dict[f] = ""

        for file_id in file_dict.keys():
            self.owner_job.fileStore.logToMaster("preparing {}".format(file_id))
            if userFileNames is None:
                temp_file_name = uuid.uuid4().hex + ".tmp"
            else:
                assert file_id in userFileNames.keys(), \
                    "Didn't find user-specified path for {}".format(file_id)
                temp_file_name = userFileNames[file_id]
            destination_path = self.work_dir + "/{fn}".format(fn=temp_file_name)
            self.owner_job.fileStore.readGlobalFile(file_id, userPath=destination_path)
            assert(os.path.exists(destination_path)), "[]"
            file_dict[file_id] = (destination_path, temp_file_name)
        
        return file_dict


def bwa_index_docker_call(job, bwa_fileId_map,
                          memory="10M", cores=1, disk="10M", 
                          bwa_docker_image="quay.io/ucsc_cgl/bwa"):
    # type: (toil.job.Job, string, string, 
    #        bool, 
    #        string, int, string, 
    #        string)
    """Builds all required indices (crufty files) that BWA needs and 
    imports them into the fileStore
    """
    def _run_bwa_index():
        bwa_index_parameters = ["index", dkr_reference_path]
        if DEBUG:
            job.fileStore.logToMaster("workDir: {}".format(localFiles.workDir()))
        # bwa docker call creates index files in the local working directory
        tlp.docker_call(tool=bwa_docker_image, 
                        parameters=bwa_index_parameters,
                        work_dir=localFiles.workDir())

        # import files to fileStore
        bwa_index_urls = ["file://" + localFiles.localFilePath(bwa_fileId_map["reference_fasta"]) +
                          suffix for suffix in bwa_index_file_suffixes()]
        new_ids        = [job.fileStore.importFile(x) for x in bwa_index_urls]
        
        # make a map of suffix files to their file Id
        index_fileId_map = dict([(suf, fid) for suf, fid in izip(bwa_index_file_suffixes(), new_ids)])

        # remove the local files
        removed = localFiles.safeLocalDelete(bwa_index_urls)
        
        if DEBUG:
            for i, x in enumerate(removed):
                if x is False:
                    job.fileStore.logToMaster("[bwa_docker_call::_run_bwa_index] " 
                                              "FAILED to remove {}".format(bwa_index_urls[i]))
                else:
                    job.fileStore.logToMaster("[bwa_docker_call::_run_bwa_index] "
                                              "removed {}".format(bwa_index_urls[i]))
        
        return index_fileId_map

    # get a local copy of the reference and reads files for docker
    localFiles = LocalFileManager(job, [bwa_fileId_map["reference_fasta"]])
    
    # arguments for the bwa indexing and alignment
    dkr_reference_path = DOCKER_DIR + localFiles.localFileName(bwa_fileId_map["reference_fasta"])
    return _run_bwa_index()


def bwa_docker_align(job, bwa_input_map, bwa_index_map, bwa_docker_image="quay.io/ucsc_cgl/bwa"):
    if DEBUG:
        job.fileStore.logToMaster("BWA input {}".format(bwa_input_map.__str__())) 
        job.fileStore.logToMaster("BWA index {}".format(bwa_index_map.__str__())) 
    
    # move the read and reference files to local working directory
    uid1, uid2, uid3 = uuid.uuid4().hex, uuid.uuid4().hex, uuid.uuid4().hex
    user_paths = {  # map of fileIds to desired file names
        bwa_input_map["reference_fasta"]   : "ref{}.fa".format(uid1),
        bwa_input_map["reads_master_fasta"]: "{}.fa".format(uid2),
    }
    for suffix in bwa_index_file_suffixes():
        user_paths[bwa_index_map[suffix]] = "ref{uid}.fa{suff}".format(uid=uid1, suff=suffix)  

    localFiles = LocalFileManager(job=job, 
                                  fileIds_to_get=bwa_input_map.values() + bwa_index_map.values(), 
                                  userFileNames=user_paths)

    job.fileStore.logToMaster("workDir %s " % localFiles.workDir())
    job.fileStore.logToMaster("ref %s " % localFiles.localFileName(bwa_input_map["reference_fasta"]))
    job.fileStore.logToMaster("sa %s " % localFiles.localFileName(bwa_index_map[".sa"]))

    dkr_reference_path = DOCKER_DIR + localFiles.localFileName(bwa_input_map["reference_fasta"])
    dkr_reads_path     = DOCKER_DIR + localFiles.localFileName(bwa_input_map["reads_master_fasta"])

    job.fileStore.logToMaster("docker ref path %s" % dkr_reference_path)
    job.fileStore.logToMaster("docker reads path %s" % dkr_reads_path)
    
    bwa_mem_parameters = ["mem", dkr_reference_path, dkr_reads_path]
    
    bwa_output_map = {}
    output_path    = localFiles.workDir() + "aln{}.sam".format(uid3)
    with open(output_path, 'w') as out_aln:
        tlp.docker_call(tool=bwa_docker_image, 
                        parameters=bwa_mem_parameters, 
                        work_dir=localFiles.workDir(), 
                        outfile=out_aln)
        assert os.path.exists(output_path)
        bwa_output_map["alignment"] = job.fileStore.importFile("file://" + output_path)
        job.fileStore.logToMaster("imported %s to %s" % (output_path, bwa_output_map["alignment"]))
        localFiles.safeLocalDelete([output_path])

    return bwa_output_map
    

def bwa_export_alignment(job, bwa_output_map, out_sam_path):
    job.fileStore.logToMaster("->>>exporting %s to %s " % (bwa_output_map["alignment"], out_sam_path))
    job.fileStore.exportFile(bwa_output_map["alignment"], out_sam_path)
    return


def bwa_docker_alignment_root(job, reference_file_id, reads_file_id, out_sam_path,
                              memory="10M", cores=1, disk="10M", 
                              bwa_docker_image="quay.io/ucsc_cgl/bwa"):
    # maps the various files needed to their unique fileStoreId, used
    # throughout the alignment pipeline
    bwa_input_map = {
        "reference_fasta": reference_file_id,
        "reads_master_fasta": reads_file_id,
    }
    bwa_index_map = job.addChildJobFn(bwa_index_docker_call, bwa_input_map).rv()
    alignment_job = job.addFollowOnJobFn(bwa_docker_align, bwa_input_map, bwa_index_map)
    bwa_output_map = alignment_job.rv()
    alignment_job.addFollowOnJobFn(bwa_export_alignment, bwa_output_map, out_sam_path)
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
