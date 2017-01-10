#!/usr/bin/env python
"""Tests for toil-nanopore
"""
from __future__ import print_function
import sys
import os
import subprocess
import textwrap
import unittest
import pysam
import numpy as np
from itertools import izip
from argparse import ArgumentParser
from margin.marginCallerLib import vcfRead
from margin.utils import ReadAlignmentStats
from margin.toil.hmm import Hmm


def baseDirectory():
    return os.getcwd()


class SubprogramCiTests(unittest.TestCase):
    def setUp(self):
        self.workdir    = os.path.join(baseDirectory(), "tests")
        self.is_file    = os.path.exists
        self.test_reads = os.path.join(self.workdir, "reads.fq")
        self.test_bam   = os.path.join(self.workdir, "bwa_bam.bam")
        self.references = os.path.join(self.workdir, "references.fa")
        self.jobstore   = "testJobStore"
        self.test_model = os.path.join(self.workdir, "last_hmm_20.txt")
        self.out_model  = os.path.join(self.workdir, "testHmm.hmm")
        self.out_sam    = os.path.join(self.workdir, "testalignment.sam")
        self.out_vcf    = os.path.join(self.workdir, "testvcf.vcf")
        self.out_stats  = os.path.join(self.workdir, "stats.txt")
        self.mut_sam    = self.getFile("inputBigMutationsBwa.bam")
        self.mut_ref    = self.getFile("referencesMutated.fa")
        self.mutations  = self.getFile("mutations.txt")

        self.assertTrue(self.is_file(self.test_reads), "Can't find reads here %s" % self.test_reads)
        self.assertTrue(self.is_file(self.references), "Can't find reference here %s" % self.references)
        self.assertTrue(self.is_file(self.test_model), "Can't find input model here %s" % self.test_model)
        #self.assertTrue(self.is_file(self.toilscript), "Can't find toil script here %s" % self.test_model)
        self.assertTrue(not self.is_file(self.out_model))
        self.assertTrue(not self.is_file(self.out_sam))
        self.assertTrue(not self.is_file(self.out_vcf))

        self.manifest_path   = None
        self.manifest_string = None
        self.config_path     = None
        self.crufty_files    = []
        self.config_string   = textwrap.dedent("""
                align:{align}
                caller:{caller}
                stats:{stats}
                debug: True
                ref:{ref}
                output_sam_path:{out_sam}
                no_chain:{no_chain}
                no_realign:{no_realign}
                hmm_file:{hmm_file}
                EM:{EM}
                output_model:{out_model}
                train_emissions:{train_emissions}
                output_vcf_path:{out_vcf}
                error_model:{error_model}
                no_margin:{no_margin}
                stats_outfile_url:{stats_out}
                set_Jukes_Cantor_emissions:
                learn_model:
                signal:
                gap_gamma:                   0.5
                match_gamma:                 0.0
                max_length_per_job:          700000
                em_iterations:               3
                model_type:                  fiveState
                max_sample_alignment_length: 50000
                random_start:                False
                update_band:                 False
                gc_content:                  0.5
                variant_threshold:           0.3
                local_alignment:             False
                noStats:                     False
                printValuePerReadAlignment:  True
                identity:                    True
                readCoverage:                True
                mismatchesPerAlignedBase:    True
                deletionsPerReadBase:        True
                insertionsPerReadBase:       True
                readLength:                  True
                """[1:])

    def getFile(self, filename):
        filepath = os.path.join(self.workdir, filename)
        self.assertTrue(self.is_file(filepath))
        return filepath

    @staticmethod
    def fileUrlify(path_string):
        return "file://" + path_string

    def tearDown(self):
        def remove_if_there(f):
            if self.is_file(f):
                os.remove(f)
        map(remove_if_there, self.crufty_files)

    def generate_manifest(self):
        self.assertTrue(self.manifest_path is not None)
        self.assertTrue(self.manifest_string is not None)

        with open(self.manifest_path, "w") as fH:
            fH.write(self.manifest_string)

        self.assertTrue(self.is_file(self.manifest_path))
        self.crufty_files.append(self.manifest_path)

    def generate_config(self, config_args):
        self.assertTrue(self.config_path is not None)
        self.assertTrue(self.config_string is not None)

        with open(self.config_path, "w") as fH:
            fH.write(self.config_string.format(**config_args))

        self.assertTrue(self.is_file(self.config_path))
        self.crufty_files.append(self.config_path)

    @staticmethod
    def base_command(command_args):
        cmd = "toil-nanopore run file:{jobstore} --config={config} --manifest={manifest} "\
              "--workDir={workdir} --clean=always".format(**command_args)
        return cmd

    def run_pipe(self):
        command_args = {
            "jobstore" : self.jobstore,
            "config"   : self.config_path,
            "manifest" : self.manifest_path,
            "workdir"  : self.workdir,
        }
        try:
            command = self.base_command(command_args)
            subprocess.check_call(command.split())
        except subprocess.CalledProcessError:
            self.assertTrue(False, "Command {} failed".format(command))

    def checkSamfile(self, observed, expected):
        """compares one sam alignment against another fails if the alignments aren't the same
        """
        self.assertTrue(self.is_file(expected))
        self.assertTrue(self.is_file(observed))

        test_sam = pysam.Samfile(observed, 'r')
        corr_sam = pysam.Samfile(expected, 'r')

        for obs, exp in izip(test_sam, corr_sam):
            self.assertTrue(obs.compare(exp) == 0)

    def validateSamfile(self, global_alignment=True):
        self.assertTrue(self.is_file(self.out_sam))
        self.assertTrue(self.is_file(self.test_reads))
        self.assertTrue(self.is_file(self.references))
        return ReadAlignmentStats.getReadAlignmentStats(samFile=self.out_sam,
                                                        readFastqFile=self.test_reads,
                                                        referenceFastaFile=self.references,
                                                        globalAlignment=global_alignment)

    def checkHmm(self):
        self.assertTrue(self.is_file(self.out_model))
        try:
            Hmm.loadHmm(self.out_model)
        except AssertionError:
            self.assertTrue(False)

    def validateVcf(self):
        mutations = set(map(lambda x : (x[0], int(x[1]) + 1, x[2]),
                            map(lambda x : x.split(), open(self.mutations, 'r'))))
        imputedMutations = vcfRead(self.out_vcf)
        intersectionSize = float(len(mutations.intersection(imputedMutations)))
        return intersectionSize / len(imputedMutations) if len(imputedMutations) else 0.0, \
            intersectionSize / len(mutations) if len(mutations) else 0.0, \
            len(imputedMutations), len(mutations)

    @staticmethod
    def print_stats(test_name, stats):
        identity   = np.mean(map(lambda x: x.identity(), stats))
        mismatches = np.mean(map(lambda x: x.mismatchesPerAlignedBase(), stats))
        insertions = np.mean(map(lambda x: x.insertionsPerReadBase(), stats))
        deletions  = np.mean(map(lambda x: x.deletionsPerReadBase(), stats))
        print("\n{test}:\n\tIdentity: {id}\n\tMismatches per aligned base: {mismat}\n\t"
              "Insertions per aligned base: {insert}\n\tDeletions per aligned base: {dels}"
              "".format(test=test_name, id=identity, mismat=mismatches, insert=insertions,
                        dels=deletions), file=sys.stdout)
        return

    def testMarginAlign(self, input_file, input_type, expected_alignment,
                        no_realign="", no_chain="", test_label=""):
        self.manifest_path   = os.path.join(self.workdir, "manifest_{}.tsv".format(test_label))
        self.manifest_string = "\t".join([input_type, (self.fileUrlify(input_file))])
        self.generate_manifest()
        config_args = {
            "align"           : " True",
            "caller"          : "",
            "stats"           : "",
            "ref"             : " {}".format(self.fileUrlify(self.references)),
            "out_sam"         : " {}".format(self.fileUrlify(self.out_sam)),
            "no_chain"        : no_chain,
            "no_realign"      : no_realign,
            "hmm_file"        : " {}".format(self.fileUrlify(self.test_model)),
            "EM"              : "",
            "out_model"       : "",
            "train_emissions" : "",
            "out_vcf"         : "",
            "error_model"     : "",
            "no_margin"       : "",
            "stats_out"       : "",
        }
        self.config_path = os.path.join(self.workdir, "config_{}.yaml".format(test_label))
        self.generate_config(config_args)
        self.crufty_files.append(self.out_sam)
        self.run_pipe()
        self.print_stats(test_label, self.validateSamfile())
        self.checkSamfile(self.out_sam, expected_alignment)

    def testMarginAlignWithFastqInput(self):
        expected_alignment = os.path.join(self.workdir, "bwa_chained_realign.sam")
        self.assertTrue(self.is_file(expected_alignment))
        self.testMarginAlign(self.test_reads, "fq", expected_alignment, test_label="testMarginAlignFastq")

    def testMarginAlignWithFastqInputNoRealign(self):
        expected_alignment = os.path.join(self.workdir, "bwa_chained.sam")
        self.assertTrue(self.is_file(expected_alignment))
        self.testMarginAlign(self.test_reads, "fq", expected_alignment, no_realign=" True",
                             test_label="testMarginAlignNoRealign")

    def testMarginAlignWithFastqInputNoChain(self):
        expected_alignment = os.path.join(self.workdir, "bwa_nochain_realign.sam")
        self.assertTrue(self.is_file(expected_alignment))
        self.testMarginAlign(self.test_reads, "fq", expected_alignment, no_chain=" True",
                             test_label="testMarginAlignNoChain")

    def testMarginAlignWithBamInput(self):
        expected_alignment = os.path.join(self.workdir, "bwa_chained_realign.sam")
        self.assertTrue(self.is_file(expected_alignment))
        self.testMarginAlign(self.test_bam, "bam", expected_alignment, test_label="testMarginAlignBam")

    def testMarginAlignWithBamInputNoRealign(self):
        expected_alignment = os.path.join(self.workdir, "bwa_chained.sam")
        self.assertTrue(self.is_file(expected_alignment))
        self.testMarginAlign(self.test_bam, "bam", expected_alignment, no_realign=" True",
                             test_label="testMarginAlignBamNoRealign")

    def testMarginAlignWithBamInputNoChain(self):
        expected_alignment = os.path.join(self.workdir, "bwa_nochain_realign.sam")
        self.assertTrue(self.is_file(expected_alignment))
        self.testMarginAlign(self.test_bam, "bam", expected_alignment, no_chain=" True",
                             test_label="testMarginAlignBamNoChain")

    def testMarginAlignEm(self, input_file, input_type, input_model,
                          no_realign="", no_chain="", test_label=""):
        self.manifest_path   = os.path.join(self.workdir, "manifest_{}.tsv".format(test_label))
        self.manifest_string = "\t".join([input_type, (self.fileUrlify(input_file))])
        self.generate_manifest()
        config_args = {
            "align"           : " True",
            "caller"          : "",
            "stats"           : "",
            "ref"             : " {}".format(self.fileUrlify(self.references)),
            "out_sam"         : " {}".format(self.fileUrlify(self.out_sam)),
            "no_chain"        : no_chain,
            "no_realign"      : no_realign,
            "hmm_file"        : input_model,
            "EM"              : " True",
            "out_model"       : " {}".format(self.fileUrlify(self.out_model)),
            "train_emissions" : " True",
            "out_vcf"         : "",
            "error_model"     : "",
            "no_margin"       : "",
            "stats_out"       : "",
        }
        self.config_path = os.path.join(self.workdir, "config_{}.yaml".format(test_label))
        self.generate_config(config_args)
        self.crufty_files.append(self.out_sam)
        self.crufty_files.append(self.out_model)
        self.run_pipe()
        self.checkHmm()

    def testMarginAlignEmWithInputModel(self):
        self.testMarginAlignEm(self.test_bam, "bam", " {}".format(self.fileUrlify(self.test_model)))

    def testMarginAlignEmWithBlankModel(self):
        self.testMarginAlignEm(self.test_bam, "bam", "")

    def testMarginCaller(self, input_file, input_type, no_margin, test_label=""):
        self.manifest_path   = os.path.join(self.workdir, "manifest_{}.tsv".format(test_label))
        self.manifest_string = "\t".join([input_type, (self.fileUrlify(input_file))])
        self.generate_manifest()
        config_args = {
            "align"           : "" ,
            "caller"          : " True",
            "stats"           : "",
            "ref"             : " {}".format(self.fileUrlify(self.mut_ref)),
            "out_sam"         : "",
            "no_chain"        : "",
            "no_realign"      : "",
            "hmm_file"        : " {}".format(self.fileUrlify(self.test_model)),
            "EM"              : "",
            "out_model"       : "",
            "train_emissions" : "",
            "out_vcf"         : " {}".format(self.fileUrlify(self.out_vcf)),
            "error_model"     : " {}".format(self.fileUrlify(self.test_model)),
            "no_margin"       : no_margin,
            "stats_out"       : "",
        }
        self.config_path = os.path.join(self.workdir, "config_{}.yaml".format(test_label))
        self.generate_config(config_args)
        self.crufty_files.append(self.out_vcf)
        self.run_pipe()

    def testMarginCallerWithMarginalization(self):
        self.testMarginCaller(self.mut_sam, "bam", "", "testMarginCaller")
        precision, recall, nCalls, nMutations = self.validateVcf()
        print("\nprecision:{precision}\nrecall:{recall}\nnCalls:{calls}\nnMutations:{muts}"
              "".format(precision=precision, recall=recall, calls=nCalls, muts=nMutations))
        return

    def testMarginCallerWithoutMarginalization(self):
        self.testMarginCaller(self.mut_sam, "bam", " True", "testMarginCallerNoMargin")
        precision, recall, nCalls, nMutations = self.validateVcf()
        print("\nprecision:{precision}\nrecall:{recall}\nnCalls:{calls}\nnMutations:{muts}"
              "".format(precision=precision, recall=recall, calls=nCalls, muts=nMutations))
        return

    def testMarginStats(self):
        self.manifest_path   = os.path.join(self.workdir, "manifest_marginStats.tsv")
        self.manifest_string = "\t".join(["bam", (self.fileUrlify(self.test_bam))])
        self.generate_manifest()
        config_args = {
            "align"           : "" ,
            "caller"          : "",
            "stats"           : " True",
            "ref"             : " {}".format(self.fileUrlify(self.references)),
            "out_sam"         : "",
            "no_chain"        : "",
            "no_realign"      : "",
            "hmm_file"        : "",
            "EM"              : "",
            "out_model"       : "",
            "train_emissions" : "",
            "out_vcf"         : "",
            "error_model"     : "",
            "no_margin"       : "",
            "stats_out"       : " {}".format(self.fileUrlify(self.out_stats)),
        }
        self.config_path = os.path.join(self.workdir, "config_marginStats.yaml")
        self.generate_config(config_args)
        self.crufty_files.append(self.out_stats)
        self.run_pipe()
        self.getFile(self.out_stats)


def main():
    testSuite = unittest.TestSuite()
    testSuite.addTest(SubprogramCiTests("testMarginAlignWithBamInput"))
    #testSuite.addTest(SubprogramCiTests("testMarginAlignWithBamInputNoRealign"))
    #testSuite.addTest(SubprogramCiTests("testMarginAlignWithBamInputNoChain"))
    #testSuite.addTest(SubprogramCiTests("testMarginAlignWithFastqInputNoRealign"))
    #testSuite.addTest(SubprogramCiTests("testMarginAlignWithFastqInputNoChain"))
    #testSuite.addTest(SubprogramCiTests("testMarginAlignWithFastqInput"))
    #testSuite.addTest(SubprogramCiTests("testMarginAlignEmWithInputModel"))
    #testSuite.addTest(SubprogramCiTests("testMarginAlignEmWithBlankModel"))
    testSuite.addTest(SubprogramCiTests("testMarginCallerWithMarginalization"))
    #testSuite.addTest(SubprogramCiTests("testMarginCallerWithoutMarginalization"))
    testSuite.addTest(SubprogramCiTests("testMarginStats"))
    
    testRunner = unittest.TextTestRunner(verbosity=1)
    testRunner.run(testSuite)

if __name__ == '__main__':
    main()
