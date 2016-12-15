import unittest
import subprocess
import os
import pysam 
from itertools import izip


class TestCase(unittest.TestCase):
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
        
        subprocess.check_call(command.split())
        
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

        subprocess.check_call(command.split())

        self.check_sam_files(test_sam_file, correct_sam_file)

        self.clean_up(test_sam_file)


def main():
    unittest.main()
        
if __name__ == '__main__':
    main()
