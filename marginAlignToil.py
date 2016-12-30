#!/usr/bin/env python
"""Toil script to run marginAlign nanopore analysis
"""
from __future__ import print_function
import os
from argparse import ArgumentParser
from toil.job import Job
from toil.common import Toil
from margin.toil.bwa import bwa_docker_alignment_root
from margin.toil.chainAlignment import chainSamFile
from margin.toil.realign import realignSamFileJobFunction
from margin.toil.expectationMaximisation import performBaumWelchOnSamJobFunction

DEBUG = False


def bwaAlignJobFunction(job, config):
    # type: (toil.job.Job, dict<string, (string and bool))>
    """Generates a SAM file, chains it (optionally), and realignes with cPecan HMM
    """
    bwa_alignment_job = job.addChildJobFn(bwa_docker_alignment_root, config)

    if bwa_alignment_job is None or config["no_chain"]:  # we're done
        return

    job.addFollowOnJobFn(chainSamFileJobFunction, config, bwa_alignment_job.rv())


def chainSamFileJobFunction(job, config, bwa_output_map):
    # Cull the files from the job store that we want
    # nb. FsId == FileStoreId
    sam_file  = job.fileStore.readGlobalFile(bwa_output_map["alignment"])
    reference = job.fileStore.readGlobalFile(config["reference_FileStoreID"])
    reads     = job.fileStore.readGlobalFile(config["sample_FileStoreID"])

    # this is the local path to the output chained SAM file
    outputSam = job.fileStore.getLocalTempFile()

    if DEBUG:
        job.fileStore.logToMaster("[chainSamFileJobFunction] samFile: {sam} output {out} "
                                  "reference {ref} reads {reads}".format(sam=sam_file,
                                                                         out=outputSam,
                                                                         ref=reference,
                                                                         reads=reads))

    chainSamFile(samFile=sam_file,
                 outputSamFile=outputSam,
                 readFastqFile=reads,
                 referenceFastaFile=reference)

    chainedSamFileId = job.fileStore.importFile("file://" + outputSam)
    job.addFollowOnJobFn(realignJobFunction, config,
                         {"chained_alignment_FileStoreID": chainedSamFileId})


def realignJobFunction(job, config, chain_alignment_output):
    if config["no_realign"]:
        job.fileStore.exportFile(chain_alignment_output["chained_alignment_FileStoreID"],
                                 config["output_sam_path"])
        return
    if config["EM"]:
        # make a child job to perform the EM and generate and import the new model
        job.fileStore.logToMaster("[realignJobFunction]Going on to run EM training "
                                  "with SAM file {sam}, read fastq {reads} and reference "
                                  "{reference}".format(
                                      sam=chain_alignment_output["chained_alignment_FileStoreID"],
                                      reads=config["sample_label"],
                                      reference=config["reference_label"]))
        assert(config["output_model"] is not None), "[realignJobFunction]Need a place to put trained_model"
        job.addChildJobFn(performBaumWelchOnSamJobFunction, config, chain_alignment_output)
        if DEBUG:
            job.fileStore.logToMaster("[realignJobFunction]Finished EM using the trained model for "
                                      "realignment")

    job.fileStore.logToMaster("[realignJobFunction]Going on to HMM realignment")
    job.addFollowOnJobFn(realignSamFileJobFunction, config, chain_alignment_output)


def main():
    def parse_args():
        parser = ArgumentParser()
        # basic input
        parser.add_argument("--reference", "-r", dest="reference", required=True)
        parser.add_argument("--reads", "-q", dest="reads", required=True)
        parser.add_argument("--hmm", dest="hmm_file", help="Hmm model parameters", required=False, default=None)
        parser.add_argument("--out_sam", "-o", dest="out_sam", required=True)
        # options
        parser.add_argument("--no_realign", dest="no_realign", help="Don't run any realignment step",
                            default=False, action="store_true")
        parser.add_argument("--no_chain", dest="no_chain", help="Don't run any chaining step",
                            default=False, action="store_true")
        parser.add_argument("--gapGamma", dest="gapGamma", help="Set the gap gamma for the AMAP function",
                            default=0.5, type=float)
        parser.add_argument("--matchGamma", dest="matchGamma", help="Set the match gamma for the AMAP function",
                            default=0.0, type=float)
        # EM options
        parser.add_argument("--em", dest="em", help="Run expectation maximisation (EM)",
                            default=False, action="store_true")
        parser.add_argument("--model_type", dest="model_type", action="store", default=None)
        parser.add_argument("--out_model", dest="out_model", action="store", default=None)
        Job.Runner.addToilOptions(parser)
        return parser.parse_args()

    args = parse_args()

    # TODO need a more streamlined way to import files
    ref_import_str    = "file://{abs_path}".format(abs_path=os.path.abspath(args.reference))
    query_import_str  = "file://{abs_path}".format(abs_path=os.path.abspath(args.reads))
    output_export_str = "file://{abs_path}".format(abs_path=os.path.abspath(args.out_sam))

    if args.em:
        assert(args.out_model is not None)
        model_export_str = "file://{abs_path}".format(abs_path=os.path.abspath(args.out_model))
    else:
        model_export_str = None

    if args.hmm_file is not None:
        hmm_import_str = "file://{abs_path}".format(abs_path=os.path.abspath(args.hmm_file))
    else:
        hmm_import_str = None

    with Toil(args) as toil:
        if not toil.options.restart:
            reference_file_id = toil.importFile(ref_import_str)
            query_file_id     = toil.importFile(query_import_str)
            hmm_file_id       = toil.importFile(hmm_import_str) if hmm_import_str is not None else None

            CONFIG = {
                # TODO consider making input model strings for logging
                "reference_label"                      : ref_import_str,
                "sample_label"                         : query_import_str,
                "reference_FileStoreID"                : reference_file_id,
                "sample_FileStoreID"                   : query_file_id,
                "input_hmm_FileStoreID"                : hmm_file_id,
                # alignment options
                "no_chain"                             : args.no_chain,
                "no_realign"                           : args.no_realign,
                "output_sam_path"                      : output_export_str,
                "gap_gamma"                            : args.gapGamma,
                "match_gamma"                          : args.matchGamma,
                # EM options
                "EM"                                   : args.em,
                "em_iterations"                        : 5,
                "random_start"                         : False,
                "set_Jukes_Cantor_emissions"           : None,
                "model_type"                           : args.model_type,
                "train_emissions"                      : True,
                "update_band"                          : False,
                "output_model"                         : model_export_str,
                "gc_content"                           : 0.5,
                "normalized_trained_model_FileStoreID" : None,
                # job-related options (sharding)
                #"max_length_per_job"                   : 5000,
                "max_length_per_job"                   : 700000,
                "max_sample_alignment_length"          : 50000,
            }

            root_job = Job.wrapJobFn(bwaAlignJobFunction, CONFIG)
            return toil.start(root_job)
        else:
            toil.restart()

if __name__ == "__main__":
    main()
