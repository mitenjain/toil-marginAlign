"""JobWrappingJobFunctions for marginCaller
"""
from __future__ import print_function
from toil_lib import require
from margin.toil.realign import shardSamJobFunction
from margin.toil.variantCaller import \
    calculateAlignedPairsJobFunction, callVariantsWithAlignedPairsJobFunction


def marginCallerJobFunction(job, config, input_samfile_fid):
    require(input_samfile_fid is not None, "[marginCallerJobFunction]input_samfile_fid is NONE")
    require(config["output_vcf_path"] is not None, "[marginCallerJobFunction]need to specify output VCF path")
    job.addFollowOnJobFn(shardSamJobFunction, config, input_samfile_fid, calculateAlignedPairsJobFunction,
                         callVariantsWithAlignedPairsJobFunction)
    # TODO add logic for folllOn stats job
