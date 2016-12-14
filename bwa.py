#!/usr/bin/env python
"""Burrow-Wheeler Aligner for pairwise alignment between DNA sequences
Module for running a BWA alignment program in a docker container.
Cite: Li H. (2013) Aligning sequence reads, clone sequences and assembly
contigs with BWA-MEM. arXiv:1303.3997v2
"""
from __future__ import print_function
import os
import uuid
from itertools import izip
import toil_lib.programs as tlp
from localFileManager import LocalFileManager

DEBUG = True
DOCKER_DIR = "/data/"


def bwa_index_file_suffixes():
    return [".amb", ".ann", ".bwt", ".pac", ".sa"]


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
