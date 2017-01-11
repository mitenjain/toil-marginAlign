"""JobWrappingJobFunctions for marginCaller
"""
from __future__ import print_function
from toil_lib import require
from margin.toil.realign import shardSamJobFunction
from margin.toil.variantCaller import \
    calculateAlignedPairsJobFunction, callVariantsWithAlignedPairsJobFunction
from margin.toil.stats import marginStatsJobFunction


def marginCallerJobFunction(job, config, input_samfile_fid, output_label):
    require(input_samfile_fid is not None, "[marginCallerJobFunction]input_samfile_fid is NONE")
    job.addFollowOnJobFn(shardSamJobFunction, config, input_samfile_fid, output_label,
                         calculateAlignedPairsJobFunction, callVariantsWithAlignedPairsJobFunction)
    if config["stats"]:
        job.addFollowOnJobFn(marginStatsJobFunction, config, input_samfile_fid, output_label)
