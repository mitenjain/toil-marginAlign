#!/usr/bin/env python2.7
from __future__ import print_function

import sys
import argparse
import os
import textwrap
import yaml
import uuid
from urlparse import urlparse

from toil.common import Toil
from toil.job import Job
from toil_lib import UserError, require
from toil_lib.files import generate_file
from toil_lib.urls import download_url_job, download_url
from toil_lib.programs import docker_call
from margin.toil.localFileManager import LocalFile
from margin.toil.alignment import AlignmentStruct, AlignmentFormat

from sample import Sample
from marginAlignToil import bwaAlignJobFunction, chainSamFileJobFunction
from marginCallerToil import marginCallerJobFunction


def getFastqFromBam(job, bam_sample, samtools_image="quay.io/ucsc_cgl/samtools"):
    # n.b. this is NOT a jobFunctionWrappingJob, it just takes the parent job as 
    # an argument to have access to the job store
    # download the BAM to the local directory, use a uid to aviod conflicts
    uid           = uuid.uuid4().hex
    work_dir      = job.fileStore.getLocalTempDir()
    local_bam     = LocalFile(workdir=work_dir, filename="bam_{}.bam".format(uid))
    bam_file_path = download_url(url=bam_sample.URL, work_dir=work_dir, name=local_bam.filenameGetter())
    fastq_reads   = LocalFile(workdir=work_dir, filename="fastq_reads{}.fq".format(uid))

    require(os.path.exists(bam_file_path), "[getFastqFromBam]didn't download BAM")
    require(bam_file_path == local_bam.fullpathGetter(), "[getFastqFromBam]bam_file_path != local_bam path")
    require(not os.path.exists(fastq_reads.fullpathGetter()), "[getFastqFromBam]fastq file already exists")

    # run samtools to get the reads from the BAM
    samtools_parameters = ["fastq", "/data/{}".format(local_bam.filenameGetter())]
    with open(fastq_reads.fullpathGetter(), 'w') as fH:
        docker_call(tool=samtools_image, parameters=samtools_parameters, work_dir=work_dir, outfile=fH)

    require(os.path.exists(fastq_reads.fullpathGetter()), "[getFastqFromBam]didn't generate reads")

    os.remove(bam_file_path)
    require(not os.path.exists(bam_file_path), "[getFastqFromBam]didn't remove bam file")

    # upload fastq to fileStore
    return job.fileStore.writeGlobalFile(fastq_reads.fullpathGetter())


def run_tool(job, config, sample):
    def cull_sample_files():
        if sample.file_type == "fq":
            config["sample_FileStoreID"] = job.addChildJobFn(download_url_job, sample.URL, disk="1G").rv()
            return None
        elif sample.file_type == "bam":
            bwa_alignment_fid = job.addChildJobFn(download_url_job, sample.URL, disk="1G").rv()
            config["sample_FileStoreID"] = job.addChildJobFn(getFastqFromBam, sample).rv()
            return bwa_alignment_fid
        else:
            require(False, "[cull_sample_files]Unsupported file type {}".format(sample.file_type))

    # download the reference
    config["reference_FileStoreID"] = job.addChildJobFn(download_url_job, config["ref"], disk="1G").rv()

    # cull the sample, which can be a fastq or a BAM
    bwa_alignment_fid = cull_sample_files()

    # download the input model, if given. Fail if no model is given and we're performing HMM realignment.
    if config["hmm_file"] is not None:
        config["input_hmm_FileStoreID"] = job.addChildJobFn(download_url_job, config["hmm_file"], disk="10M").rv()
    else:
        if not config["no_realign"]:
            require(config["EM"], "[run_tool]Need to specify an input model or set EM to "
                                  "True to perform HMM realignment")
        config["input_hmm_FileStoreID"] = None

    # initialize key in config for trained model if we're performing EM
    if config["EM"]:
        config["normalized_trained_model_FileStoreID"] = None

    # TODO refactor this out
    config["sample_label"]    = sample.URL
    config["reference_label"] = config["ref"]

    # Pipeline starts here
    if config["align"]:
        if bwa_alignment_fid is None:
            marginAlignJob = job.addFollowOnJobFn(bwaAlignJobFunction, config)
        else:
            aln_struct = AlignmentStruct(bwa_alignment_fid, AlignmentFormat.BAM)
            marginAlignJob = job.addFollowOnJobFn(chainSamFileJobFunction, config, aln_struct)
        # TODO add followon to marginCaller here.
    if config["learn_model"]:
        raise NotImplementedError
    if config["caller"]:
        require(config["error_model"] is not None, "[run_tool]Caller subprogram requires an error model"
                "None was given")
        config["error_model_FileStoreID"] = job.addChildJobFn(download_url_job,
                                                              config["error_model"],
                                                              disk="10M").rv()
        job.addFollowOnJobFn(marginCallerJobFunction, config, bwa_alignment_fid)
    if config["stats"]:
        raise NotImplementedError
    if config["signal"]:
        raise NotImplementedError


def print_help():
    """this is the help
    """
    return print_help.__doc__


def generateConfig():
    return textwrap.dedent("""
        # UCSC Nanopore Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline: "toil-nanopore run"
        #
        # URLs can take the form: http://, ftp://, file://, s3://, gnos://
        # Local inputs follow the URL convention: file:///full/path/to/input
        # S3 URLs follow the convention: s3://bucket/directory/file.txt
        #
        # some options have been filled in with defaults

        ## Universal Options/Inputs ##
        # Required: Which subprograms to run, typically you run all 4, but you can run them piecemeal if you like
        # in that case the provided inputs will be checked at run time
        align:
        learn_model:
        caller: True
        stats:
        signal:

        # Optional: Debug increasing logging
        debug: True

        # Required: Reference fasta file
        ref: file:///Users/Rand/projects/toil_dev/toil-marginAlign/tests/references.fa

        ## MarginAlign Options ##
        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        output_sam_path:

        # all required options have default values
        gap_gamma:                     0.5
        match_gamma:                   0.0
        "max_length_per_job":          700000
        "max_sample_alignment_length": 50000

        # Optional: no chain and/or no re-align
        no_chain:
        no_realign:

        # Optional: Alignment Model
        hmm_file: file:///Users/Rand/projects/toil_dev/toil-marginAlign/tests/last_hmm_20.txt

        # EM options
        # Optional: set true to do EM
        EM:
        # Required: path to location of output model
        output_model:
        model_type:    fiveState
        em_iterations: 5
        random_start:  False
        # set_Jukes_Cantor_emissions is of type Float
        set_Jukes_Cantor_emissions:
        update_band:     False
        gc_content:      0.5
        train_emissions: True

        ## MarginCaller Options ##
        # Required: path to VCF output
        output_vcf_path: file:///Users/Rand/projects/toil_dev/toil-marginAlign/sandbox/testvcf.vcf
        # Required: Error model
        error_model: file:///Users/Rand/projects/toil_dev/toil-marginAlign/tests/last_hmm_20.txt

        # Options
        # n.b. required options have default values filled in
        no_margin: False
        variant_threshold: 0.3
    """[1:])


def generateManifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information for each sample to be run.
        #
        #   Lines should contain two tab-seperated fields: file_type and URL
        #   file_type options:
        #       fq-gzp gzipped file of read sequences in FASTQ format
        #           fq file of read sequences in FASTQ format
        #          bam alignment file in BAM format (sorted or unsorted)
        #       fa-gzp gzipped file of read sequences in FASTA format
        #           fa file of read sequences in FASTA format
        #       f5-tar tarball of MinION, basecalled, .fast5 files
        #   NOTE: as of 1/3/16 only fq implemented
        #   Eg:
        #   fq-tar  file://path/to/file/reads.tar
        #   f5-tar  s3://my-bucket/directory/
        #   bam     file://path/to/giantbam.bam
        #   Place your samples below, one sample per line.
        """[1:])


def parseManifest(path_to_manifest):
    require(os.path.exists(path_to_manifest), "[parseManifest]Didn't find manifest file, looked "
            "{}".format(path_to_manifest))
    allowed_file_types = ["fq-gzp", "fq", "fa-gzp", "fa", "f5-tar", "bam"]

    def parse_line(line):
        # double check input, shouldn't need to though
        require(not line.isspace() and not line.startswith("#"), "[parse_line]Invalid {}".format(line))
        sample = line.strip().split("\t")
        # there should only be two entries, the file_type and the URL
        require(len(sample) == 2, "[parse_line]Invalid, len(line) != 2, offending {}".format(line))
        file_type, sample_url = sample
        # check the file_type and the URL
        require(file_type in allowed_file_types, "[parse_line]Unrecognized file type {}".format(file_type))
        require(urlparse(sample_url).scheme and urlparse(sample_url), "Invalid URL passed for {}".format(sample_url))
        return Sample(file_type=file_type, URL=sample_url, file_id=None)

    with open(path_to_manifest, "r") as fH:
        return map(parse_line, [x for x in fH if (not x.isspace() and not x.startswith("#"))])


def main():
    """toil-nanopore master script
    """
    def parse_args():
        parser = argparse.ArgumentParser(description=print_help.__doc__,
                                         formatter_class=argparse.RawTextHelpFormatter)
        subparsers = parser.add_subparsers(dest="command")
        run_parser = subparsers.add_parser("run", help="runs nanopore pipeline with config")
        subparsers.add_parser("generate", help="generates a config file for your run, do this first")

        group = run_parser.add_mutually_exclusive_group()
        group.add_argument('--config', default='config-toil-nanopore.yaml', type=str,
                           help='Path to the (filled in) config file, generated with "generate".')
        group.add_argument('--manifest', default='manifest-toil-nanopore.tsv', type=str,
                           help='Path to the (filled in) manifest file, generated with "generate". '
                                '\nDefault value: "%(default)s".')
        Job.Runner.addToilOptions(run_parser)

        return parser.parse_args()

    def exitBadInput(message=None):
        if message is not None:
            print(message, file=sys.stderr)
        sys.exit(1)

    if len(sys.argv) == 1:
        exitBadInput(print_help())

    cwd = os.path.dirname(os.path.abspath(__file__))

    args = parse_args()

    if args.command == "generate":
        try:
            config_path = os.path.join(cwd, "config-toil-nanopore.yaml")
            generate_file(config_path, generateConfig)
        except UserError:
            print("[toil-nanopore]NOTICE using existing config file {}".format(config_path))
            pass
        try:
            manifest_path = os.path.join(cwd, "manifest-toil-nanopore.tsv")
            generate_file(manifest_path, generateManifest)
        except UserError:
            print("[toil-nanopore]NOTICE using existing manifest {}".format(manifest_path))

    elif args.command == "run":
        require(os.path.exists(args.config), "{config} not found run generate-config".format(config=args.config))
        # Parse config
        config  = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        samples = parseManifest(args.manifest)
        for sample in samples:
            with Toil(args) as toil:
                if not toil.options.restart:
                    root_job = Job.wrapJobFn(run_tool, config, sample)
                    return toil.start(root_job)
                else:
                    toil.restart()


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
