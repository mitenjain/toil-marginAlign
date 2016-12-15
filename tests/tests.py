import unittest
import subprocess
import os


class TestCase(unittest.TestCase):
    def test_bwa(self):
        command = "python marginAlignToil.py file:jobstore "\
                  "-r ./tests/references.fa -q ./tests/reads.fq "\
                  "-o ./testchained.sam --workDir={cwd} --logInfo --no_chain".format(cwd=os.getcwd())
        subprocess.check_call(command.split())


def main():
    unittest.main()
        
if __name__ == '__main__':
    main()
