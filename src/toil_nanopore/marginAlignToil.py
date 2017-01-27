#!/usr/bin/env python
"""Toil script to run marginAlign analysis
"""
from __future__ import print_function
import os
from argparse import ArgumentParser
from toil.job import Job
from toil.common import Toil

from margin.toil.bwa import bwa_docker_alignment_root
from margin.toil.realign import realignSamFileJobFunction
from margin.toil.chainAlignment import chainSamFile
from margin.toil.localFileManager import LocalFile, deliverOutput
from margin.toil.expectationMaximisation import performBaumWelchOnSamJobFunction


def baseDirectoryPath():
    return os.path.dirname(os.path.abspath(__file__)) + "/"


def getDefaultHmm():
    default_model = baseDirectoryPath() + "tests/last_hmm_20.txt"
    assert(os.path.exists(default_model)), "[getDefaultHmm]ERROR didn't find default model, "\
                                           "looked for it {here}".format(here=default_model)
    return default_model


def bwaAlignJobFunction(job, config):
    # type: (toil.job.Job, dict<string, (string and bool))>
    """Generates a SAM file, chains it (optionally), and realignes with cPecan HMM
    """
    bwa_alignment_job = job.addChildJobFn(bwa_docker_alignment_root, config)

    job.addFollowOnJobFn(chainSamFileJobFunction, config, bwa_alignment_job.rv())


def chainSamFileJobFunction(job, config, aln_struct):
    # Cull the files from the job store that we want
    if config["chain"] is None and config["realign"] is None:
        job.fileStore.logToMaster("[chainSamFileJobFunction]Nothing to do.")
        return

    if config["chain"] is not None:
        sam_file   = job.fileStore.readGlobalFile(aln_struct.FileStoreID())
        reference  = job.fileStore.readGlobalFile(config["reference_FileStoreID"])
        reads      = job.fileStore.readGlobalFile(config["sample_FileStoreID"])
        workdir    = job.fileStore.getLocalTempDir()
        output_sam = LocalFile(workdir=workdir, filename="{}_chained.bam".format(config["sample_label"]))

        if config["debug"]:
            job.fileStore.logToMaster("[chainSamFileJobFunction] chaining {bwa_out} (locally: {sam})"
                                      "".format(bwa_out=aln_struct.FileStoreID(), sam=sam_file))

        chainSamFile(parent_job=job,
                     samFile=sam_file,
                     outputSamFile=output_sam.fullpathGetter(),
                     readFastqFile=reads,
                     referenceFastaFile=reference)

        chainedSamFileId = job.fileStore.writeGlobalFile(output_sam.fullpathGetter())
        deliverOutput(job, output_sam, config["output_dir"])
        job.addFollowOnJobFn(realignJobFunction, config, chainedSamFileId)

    else:
        job.fileStore.logToMaster("[chainSamFileJobFunction]Not chaining SAM, passing {sam} "
                                  "on to realignment".format(sam=aln_struct.FileStoreID()))
        job.addFollowOnJobFn(realignJobFunction, config, aln_struct.FileStoreID())


def realignJobFunction(job, config, input_samfile_fid):
    if config["realign"] is None:  # the chained SAM has already been delivered
        return
    if config["EM"]:
        # make a child job to perform the EM and generate and import the new model
        job.fileStore.logToMaster("[realignJobFunction]Queueing EM training "
                                  "with SAM file {sam}, read fastq {reads} and reference "
                                  "{reference}".format(
                                      sam=input_samfile_fid,
                                      reads=config["sample_label"],
                                      reference=config["reference_label"]))
        job.addChildJobFn(performBaumWelchOnSamJobFunction, config, input_samfile_fid)

    # need disk = ref_size + bam_size
    # memory = ref_size(hash) + sam
    job.fileStore.logToMaster("[realignJobFunction]Queueing up HMM realignment")
    realign_label = "realigned" if config["chain"] else "noChain_realigned"
    return job.addFollowOnJobFn(realignSamFileJobFunction, config, input_samfile_fid, realign_label).rv()
