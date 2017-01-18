"""JobWrappingJobFunctions for marginCaller
"""
from __future__ import print_function
from toil_lib import require
from margin.toil.realign import shardSamJobFunction
from margin.toil.variantCaller import\
    calculateAlignedPairsJobFunction,\
    callVariantsWithAlignedPairsJobFunction1,\
    marginalizePosteriorProbsJobFunction
from margin.toil.alignment import splitLargeAlignment
from margin.toil.stats import marginStatsJobFunction


def marginCallerJobFunction(job, config, input_samfile_fid, output_label):
    require(input_samfile_fid is not None, "[marginCallerJobFunction]input_samfile_fid is NONE")
    # split up the large alignment
    smaller_alns = splitLargeAlignment(job, config, input_samfile_fid)
    expectations = []
    # this loop runs through the smaller alignments and sets a child job to get the aligned pairs for 
    # each one. then it marginalizes over the columns in the alignment and adds a promise of a dict containing
    # the posteriors to the list `expctations`
    for aln in smaller_alns:
        disk   = aln.size + config["reference_FileStoreID"].size
        memory = (6 * aln.size)
        position_expectations = job.addChildJobFn(shardSamJobFunction,
                                                  config, aln,
                                                  calculateAlignedPairsJobFunction,
                                                  marginalizePosteriorProbsJobFunction,
                                                  disk=disk, memory=memory).rv()
        expectations.append(position_expectations)

    job.addFollowOnJobFn(callVariantsWithAlignedPairsJobFunction1, config, expectations, output_label)

    if config["stats"]:
        job.addFollowOnJobFn(marginStatsJobFunction, config, input_samfile_fid, output_label,
                             memory=input_samfile_fid.size)
