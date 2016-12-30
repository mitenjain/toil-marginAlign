#!/usr/bin/env python
"""CI and unit tests for toil-marginAlign
"""
import unittest
import subprocess
import os
import pysam
from itertools import izip
from cPecan.cPecanEm import Hmm


class ToilMarginAlignCiTest(unittest.TestCase):
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
    def toil_clean(job_store):
        cmd = "toil clean {}/".format(job_store)
        subprocess.check_call(cmd.split())

    def check_sam_files(self, observed_sam, expected_sam):
        self.assertTrue(os.path.exists(expected_sam))
        self.assertTrue(os.path.exists(observed_sam))

        test_sam = pysam.Samfile(observed_sam, 'r')
        corr_sam = pysam.Samfile(expected_sam, 'r')

        for obs, exp in izip(test_sam, corr_sam):
            self.assertTrue(obs.compare(exp) == 0)

    def test_bwa(self):
        test_sam_file    = "./tests/TESTSAM.sam"
        correct_sam_file = "./tests/bwa_only.sam"

        self.initialize(test_sam_file, correct_sam_file)

        command = "python marginAlignToil.py file:jobstore "\
                  "-r ./tests/references.fa -q ./tests/reads.fq "\
                  "-o {test_sam} --workDir={cwd} --logInfo --no_chain".format(cwd=os.getcwd(),
                                                                              test_sam=test_sam_file)
        try:
            subprocess.check_call(command.split())
        except subprocess.CalledProcessError:
            self.toil_clean("jobstore")
            self.assertTrue(False)

        # test that the AlignedRegions in the observed and expected are the same
        self.check_sam_files(test_sam_file, correct_sam_file)

        self.clean_up(test_sam_file)

    def test_bwa_chained(self):
        test_sam_file    = "./tests/TESTSAM.sam"
        correct_sam_file = "./tests/bwa_chained.sam"

        self.initialize(test_sam_file, correct_sam_file)

        command = "python marginAlignToil.py file:jobstore "\
                  "-r ./tests/references.fa -q ./tests/reads.fq "\
                  "-o {test_sam} --workDir={cwd} --logInfo --no_realign".format(cwd=os.getcwd(),
                                                                                test_sam=test_sam_file)
        try:
            subprocess.check_call(command.split())
        except subprocess.CalledProcessError:
            self.toil_clean("jobstore")
            self.assertTrue(False)

        self.check_sam_files(test_sam_file, correct_sam_file)

        self.clean_up(test_sam_file)

    def test_bwa_realign(self):
        test_sam_file    = "./tests/TESTSAM.sam"
        correct_sam_file = "./tests/bwa_realigned.sam"

        self.initialize(test_sam_file, correct_sam_file)

        command = "python marginAlignToil.py file:jobstore "\
                  "-r ./tests/references.fa -q ./tests/reads.fq "\
                  "-o {test_sam} --workDir={cwd} --logInfo "\
                  "--hmm ../marginAlign/src/margin/mappers/last_hmm_20.txt".format(cwd=os.getcwd(),
                                                                                   test_sam=test_sam_file)

        try:
            subprocess.check_call(command.split())
        except subprocess.CalledProcessError:
            self.toil_clean("jobstore")
            self.assertTrue(False)

        self.check_sam_files(test_sam_file, correct_sam_file)

        self.clean_up(test_sam_file)

    def test_bwa_em(self):
        def check_hmm(hmm_file):
            self.assertTrue(os.path.exists(hmm_file))
            try:
                Hmm.loadHmm(hmm_file)
            except AssertionError:
                self.assertTrue(False)

        test_sam_file   = "./tests/TESTSAM.sam"
        test_model_file = "./tests/testmodel.hmm"

        command = "python marginAlignToil.py file:jobstore "\
                  "--workDir={cwd} -r ./tests/references.fa "\
                  "-q ./tests/reads.fq -o {test_sam} "\
                  "--em --model_type=fiveState --out_model {test_model}"\
                  "".format(cwd=os.getcwd(),
                            test_sam=test_sam_file,
                            test_model=test_model_file)

        try:
            subprocess.check_call(command.split())
        except subprocess.CalledProcessError:
            self.toil_clean("jobstore")
            self.assertTrue(False)

        check_hmm(test_model_file)
        self.clean_up(test_sam_file)
        self.clean_up(test_model_file)


def main():
    testSuite = unittest.TestSuite()
    #testSuite.addTest(ToilMarginAlignCiTest("test_bwa"))
    #testSuite.addTest(ToilMarginAlignCiTest("test_bwa_chained"))
    #testSuite.addTest(ToilMarginAlignCiTest("test_bwa_realign"))
    testSuite.addTest(ToilMarginAlignCiTest("test_bwa_em"))
    testRunner = unittest.TextTestRunner(verbosity=2)
    testRunner.run(testSuite)


if __name__ == '__main__':
    main()
