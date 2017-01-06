#!/usr/bin/env python
"""Tests for toil-nanopore
"""
from __future__ import print_function
import os
import subprocess
import textwrap
import unittest
from argparse import ArgumentParser


def baseDirectory():
    return os.getcwd()


class SubprogramCiTests(unittest.TestCase):
    def setUp(self):
        self.workdir    = os.path.join(baseDirectory(), "tests")
        self.toilscript = os.path.join(baseDirectory(), "toil-nanopore.py")
        self.test_reads = os.path.join(self.workdir, "reads.fq")
        self.test_bam   = os.path.join(self.workdir, "bwa_chained_realign.bam")
        self.references = os.path.join(self.workdir, "references.fa")
        self.jobstore   = "testJobStore"
        self.test_model = os.path.join(self.workdir, "last_hmm_20.txt")
        self.out_model  = os.path.join(self.workdir, "testHmm.hmm")
        self.out_sam    = os.path.join(self.workdir, "testalignment.sam")
        self.is_file    = os.path.exists

        self.assertTrue(self.is_file(self.test_reads), "Can't find reads here %s" % self.test_reads)
        self.assertTrue(self.is_file(self.references), "Can't find reference here %s" % self.references)
        self.assertTrue(self.is_file(self.test_model), "Can't find input model here %s" % self.test_model)
        self.assertTrue(self.is_file(self.toilscript), "Can't find toil script here %s" % self.test_model)
        self.assertTrue(not self.is_file(self.out_model))
        self.assertTrue(not self.is_file(self.out_sam))

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
                learn_model:
                signal:
                gap_gamma:                   0.5
                match_gamma:                 0.5
                max_length_per_job:          700000
                em_iterations:               5
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
        cmd = "{script} run file:{jobstore} --config={config} --manifest={manifest} "\
              "--workDir={workdir} --clean=always".format(**command_args)
        return cmd

    def run_pipe(self):
        command_args = {
            "script"   : self.toilscript,
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

    def testMarginAlign(self):
        self.manifest_path   = os.path.join(self.workdir, "manifest_testMarginAlign.tsv")
        self.manifest_string = "\t".join(["fq", ("file://" + self.test_reads)])
        self.generate_manifest()
        config_args = {
            "align"           : " True",
            "caller"          : "",
            "stats"           : "",
            "ref"             : " {}".format(self.fileUrlify(self.references)),
            "out_sam"         : " {}".format(self.fileUrlify(self.out_sam)),
            "no_chain"        : "",
            "no_realign"      : "",
            "hmm_file"        : " {}".format(self.fileUrlify(self.test_model)),
            "EM"              : "",
            "out_model"       : "",
            "train_emissions" : "",
            "out_vcf"         : "",
            "error_model"     : "",
            "no_margin"       : "",
            "stats_out"       : "",
        }
        self.config_path = os.path.join(self.workdir, "config_testMarginAlign.yaml")
        self.generate_config(config_args)
        self.run_pipe()




def main():
    parser = ArgumentParser()
    parser.add_argument("--work_dir", action="store", dest="work_dir", default="./", required=False)
    parser.add_argument("--show_stats", action="store_true", dest="show_stats", default=False, required=False)
    parser.add_argument("--debug", action="store_true", dest="debug", default=False, required=False)
    args = parser.parse_args()
    
    testSuite = unittest.TestSuite()
    testSuite.addTest(SubprogramCiTests("testMarginAlign"))
    testRunner = unittest.TextTestRunner(verbosity=1)
    testRunner.run(testSuite)

if __name__ == '__main__':
    main()
