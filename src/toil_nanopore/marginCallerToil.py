"""JobWrappingJobFunctions for marginCaller
"""
from __future__ import print_function
from toil_lib import require
from margin.toil.realign import shardSamJobFunction
from margin.toil.variantCaller import\
    calculateAlignedPairsJobFunction,\
    marginalizePosteriorProbsJobFunction
from margin.toil.variantCall import makeVcfFromVariantCallsJobFunction2
from margin.toil.alignment import shardAlignmentByRegion
from margin.toil.stats import marginStatsJobFunction
from margin.toil.hmm import downloadHmm


def marginCallerJobFunction(job, config, input_samfile_fid, output_label):
    require(input_samfile_fid is not None, "[marginCallerJobFunction]input_samfile_fid is NONE")
    # split up the large alignment by region
    smaller_alns        = shardAlignmentByRegion(job, config, input_samfile_fid)
    all_variant_calls   = []
    hidden_markov_model = downloadHmm(job, config)
    # this loop runs through the smaller alignments and sets a child job to get the aligned pairs for 
    # each one. then it marginalizes over the columns in the alignment and adds a promise of a dict containing
    # the posteriors to the list `expctations`
    job.fileStore.logToMaster("[marginCallerJobFunction]Have {} smaller alignments to variant "
                              "call".format(len(smaller_alns)))
    for aln in smaller_alns:
        disk          = input_samfile_fid.size + config["reference_FileStoreID"].size
        memory        = (6 * input_samfile_fid.size)
        variant_calls = job.addChildJobFn(shardSamJobFunction,
                                          config, aln, hidden_markov_model,
                                          calculateAlignedPairsJobFunction,
                                          marginalizePosteriorProbsJobFunction,
                                          disk=disk, memory=memory).rv()
        all_variant_calls.append(variant_calls)

    job.addFollowOnJobFn(makeVcfFromVariantCallsJobFunction2, config, all_variant_calls, output_label)

    if config["stats"]:
        job.addFollowOnJobFn(marginStatsJobFunction, config, input_samfile_fid, output_label,
                             memory=input_samfile_fid.size)
