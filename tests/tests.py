#!/usr/bin/env python
"""CI and unit tests for toil-marginAlign
"""
from __future__ import print_function
import unittest
import subprocess
import os
import sys
import pysam
import numpy as np
from itertools import izip
from argparse import ArgumentParser
from margin.toil.hmm import Hmm
from margin.utils import ReadAlignmentStats


DEVNULL = open(os.devnull, 'w')


class ToilMarginAlignCiTest(unittest.TestCase):
    def __init__(self, test_name, work_dir, show_stats, debug):
        # type: (str, str, bool, bool)
        super(ToilMarginAlignCiTest, self).__init__(test_name)
        self.work_dir   = work_dir
        self.show_stats = show_stats
        self.test_name  = test_name
        self.debug      = debug

    def setUp(self):
        self.reads_fastq     = "./tests/reads.fq"
        self.references      = "./tests/references.fa"
        self.test_jobstore   = "testjobstore"
        self.test_samfile    = "./tests/TESTSAM.sam"
        self.test_model      = "./tests/last_hmm_20.txt"
        self.test_model_file = "./tests/testmodel.hmm"
        self.assertTrue(os.path.exists(self.reads_fastq))
        self.assertTrue(os.path.exists(self.test_model))
        self.assertTrue(os.path.exists(self.references))
        self.assertTrue(not os.path.exists(self.test_samfile))
        self.assertTrue(not os.path.exists(self.test_model_file))

    @staticmethod
    def initialize(test_sam_file, correct_sam_file):
        assert(os.path.exists(correct_sam_file))
        if os.path.exists(test_sam_file):
            os.remove(test_sam_file)
        return

    @staticmethod
    def clean_up(test_sam_file):
        if os.path.exists(test_sam_file):
            os.remove(test_sam_file)
        else:
            pass
        return

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

    def toil_clean(self, job_store):
        cmd = "toil clean {}/".format(self.test_jobstore)
        subprocess.check_call(cmd.split())

    def check_sam_files(self, observed_sam, expected_sam):
        """compares one sam alignment against another fails if the alignments aren't the same
        """
        self.assertTrue(os.path.exists(expected_sam))
        self.assertTrue(os.path.exists(observed_sam))

        test_sam = pysam.Samfile(observed_sam, 'r')
        corr_sam = pysam.Samfile(expected_sam, 'r')

        for obs, exp in izip(test_sam, corr_sam):
            self.assertTrue(obs.compare(exp) == 0)

    def run_test(self, command, outmodel=None):
        try:
            subprocess.check_call(command.split(),
                                  stdout=(None if self.debug else DEVNULL),
                                  stderr=(None if self.debug else DEVNULL))
        except subprocess.CalledProcessError:
            self.toil_clean(self.test_jobstore)
            self.clean_up(self.test_samfile)
            if outmodel is not None:
                self.clean_up(outmodel)
            self.assertTrue(False)

    def validate_sam(self, samfile, reads_fastq, reference_fasta, global_alignment=True):
        self.assertTrue(os.path.exists(samfile))
        self.assertTrue(os.path.exists(reads_fastq))
        self.assertTrue(os.path.exists(reference_fasta))
        return ReadAlignmentStats.getReadAlignmentStats(samFile=samfile,
                                                        readFastqFile=reads_fastq,
                                                        referenceFastaFile=reference_fasta,
                                                        globalAlignment=global_alignment)

    def check_hmm(self):
        self.assertTrue(os.path.exists(self.test_model_file))
        try:
            Hmm.loadHmm(self.test_model_file)
        except AssertionError:
            self.assertTrue(False)

    def test_bwa_only(self):
        correct_sam_file = "./tests/bwa_only.sam"

        self.initialize(self.test_samfile, correct_sam_file)

        command = "python marginAlignToil.py file:{jobstore} "\
                  "-r {references} -q {reads} "\
                  "-o {test_sam} --workDir={cwd} --logInfo --no_chain "\
                  "--no_realign".format(jobstore=self.test_jobstore,
                                        references=self.references,
                                        reads=self.reads_fastq,
                                        cwd=os.getcwd(),
                                        test_sam=self.test_samfile)
        self.run_test(command)

        # test that the AlignedRegions in the observed and expected are the same
        self.check_sam_files(self.test_samfile, correct_sam_file)

        stats = self.validate_sam(self.test_samfile, self.reads_fastq, self.references)
        if self.show_stats:
            self.print_stats(self.test_name, stats)

        self.clean_up(self.test_samfile)

    def test_bwa_chained(self):
        correct_sam_file = "./tests/bwa_chained.sam"

        self.initialize(self.test_samfile, correct_sam_file)

        command = "python marginAlignToil.py file:{jobstore} "\
                  "-r {references} -q {reads} "\
                  "-o {test_sam} --workDir={cwd} --logInfo --no_realign".format(jobstore=self.test_jobstore,
                                                                                references=self.references,
                                                                                reads=self.reads_fastq,
                                                                                cwd=os.getcwd(),
                                                                                test_sam=self.test_samfile)

        self.run_test(command)

        self.check_sam_files(self.test_samfile, correct_sam_file)
        stats = self.validate_sam(self.test_samfile, self.reads_fastq, self.references)
        if self.show_stats:
            self.print_stats(self.test_name, stats)

        self.clean_up(self.test_samfile)

    def test_bwa_realign(self):
        correct_sam_file = "./tests/bwa_chained_realign.sam"

        self.initialize(self.test_samfile, correct_sam_file)

        command = "python marginAlignToil.py file:{jobstore} "\
                  "-r {references} -q {reads} "\
                  "-o {test_sam} --workDir={cwd} --logInfo "\
                  "--hmm {hmm} ".format(jobstore=self.test_jobstore,
                                        references=self.references,
                                        reads=self.reads_fastq,
                                        cwd=os.getcwd(),
                                        test_sam=self.test_samfile,
                                        hmm=self.test_model)

        self.run_test(command)

        self.check_sam_files(self.test_samfile, correct_sam_file)
        stats = self.validate_sam(self.test_samfile, self.reads_fastq, self.references)
        if self.show_stats:
            self.print_stats(self.test_name, stats)

        self.clean_up(self.test_samfile)

    def test_bwa_realign_no_chain(self):
        correct_sam_file = "./tests/bwa_nochain_realign.sam"

        self.initialize(self.test_samfile, correct_sam_file)

        command = "python marginAlignToil.py file:{jobstore} "\
                  "-r {references} -q {reads} "\
                  "-o {test_sam} --workDir={cwd} --logInfo --no_chain "\
                  "--hmm {hmm} ".format(jobstore=self.test_jobstore,
                                        references=self.references,
                                        reads=self.reads_fastq,
                                        cwd=os.getcwd(),
                                        test_sam=self.test_samfile,
                                        hmm=self.test_model)

        self.run_test(command)
        self.check_sam_files(self.test_samfile, correct_sam_file)
        stats = self.validate_sam(self.test_samfile, self.reads_fastq, self.references)
        if self.show_stats:
            self.print_stats(self.test_name, stats)

        self.clean_up(self.test_samfile)

    def test_defaults(self):
        raise NotImplementedError

    def test_bwa_em_no_chain(self):
        command = "python marginAlignToil.py file:{jobstore} "\
                  "--workDir={cwd} -r {reference} "\
                  "-q {reads} -o {test_sam} --no_chain "\
                  "--em --model_type=fiveState --out_model {test_model}"\
                  "".format(jobstore=self.test_jobstore,
                            reference=self.references,
                            reads=self.reads_fastq,
                            cwd=os.getcwd(),
                            test_sam=self.test_samfile,
                            test_model=self.test_model_file)

        self.run_test(command=command, outmodel=self.test_model_file)
        self.check_hmm()
        stats = self.validate_sam(self.test_samfile, self.reads_fastq, self.references)
        if self.show_stats:
            self.print_stats(self.test_name, stats)
        self.clean_up(self.test_samfile)
        self.clean_up(self.test_model_file)

    def test_bwa_em(self):
        command = "python marginAlignToil.py file:{jobstore} "\
                  "--workDir={cwd} -r {reference} "\
                  "-q {reads} -o {test_sam} "\
                  "--em --model_type=fiveState --out_model {test_model}"\
                  "".format(jobstore=self.test_jobstore,
                            reference=self.references,
                            reads=self.reads_fastq,
                            cwd=os.getcwd(),
                            test_sam=self.test_samfile,
                            test_model=self.test_model_file)

        self.run_test(command=command, outmodel=self.test_model_file)
        self.check_hmm()
        stats = self.validate_sam(self.test_samfile, self.reads_fastq, self.references)
        if self.show_stats:
            self.print_stats(self.test_name, stats)
        self.clean_up(self.test_samfile)
        self.clean_up(self.test_model_file)


def main():
    parser = ArgumentParser()
    parser.add_argument("--work_dir", action="store", dest="work_dir", default="./", required=False)
    parser.add_argument("--show_stats", action="store_true", dest="show_stats", default=False, required=False)
    parser.add_argument("--debug", action="store_true", dest="debug", default=False, required=False)
    args = parser.parse_args()

    work_dir   = args.work_dir
    show_stats = args.show_stats
    debug      = args.debug

    testSuite = unittest.TestSuite()
    testSuite.addTest(ToilMarginAlignCiTest("test_bwa_only", work_dir, show_stats, debug))
    testSuite.addTest(ToilMarginAlignCiTest("test_bwa_chained", work_dir, show_stats, debug))
    testSuite.addTest(ToilMarginAlignCiTest("test_bwa_realign_no_chain", work_dir, show_stats, debug))
    testSuite.addTest(ToilMarginAlignCiTest("test_bwa_realign", work_dir, show_stats, debug))
    testSuite.addTest(ToilMarginAlignCiTest("test_bwa_em", work_dir, show_stats, debug))
    testSuite.addTest(ToilMarginAlignCiTest("test_bwa_em_no_chain", work_dir, show_stats, debug))
    testRunner = unittest.TextTestRunner(verbosity=2)
    testRunner.run(testSuite)


if __name__ == '__main__':
    main()
