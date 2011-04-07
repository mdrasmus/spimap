/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree evolution model

=============================================================================*/


#include "common.h"
#include "branch_prior.h"
#include "logging.h"
#include "Matrix.h"
#include "model.h"
#include "phylogeny.h"
#include "seq_likelihood.h"
#include "top_prior.h"



namespace spidir {


SpimapModel::SpimapModel(
    int nnodes, SpeciesTree *stree, 
    SpidirParams *params, 
    int *gene2species,
    float predupprob, float dupprob, float lossprob, 
    int nsamples, bool approx, bool useBranchPrior) :
    
    Model(nnodes),
    nnodes(nnodes),
    stree(stree),
    params(params),
    gene2species(gene2species),
    predupprob(predupprob),
    dupprob(dupprob),
    lossprob(lossprob),
    nsamples(nsamples),
    approx(approx),
    useBranchPrior(useBranchPrior)
{
    doomtable = new double [stree->nnodes];
    calcDoomTable(stree, dupprob, lossprob, doomtable);
}


SpimapModel::~SpimapModel()
{
    delete [] doomtable;

    if (likelihoodFunc)
        delete likelihoodFunc;
}

void SpimapModel::setTree(Tree *_tree)
{
    tree = _tree;
    spidir::reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
}


double SpimapModel::likelihood()
{
    double logp = 0.0;
    if (likelihoodFunc) {
        Timer timer;
        logp = likelihoodFunc->findLengths(tree);
        seq_runtime += timer.time();
    }
    return logp;
}

double SpimapModel::branchPrior()
{
    if (!useBranchPrior)
        return 0.0;

    Timer timer;
    const float generate = -99; // integrate over gene rate
    double logp = spidir::branchPrior(tree, stree,
                                      recon, events, params,
                                      generate, predupprob, dupprob, lossprob,
                                      nsamples, approx);
    branch_runtime += timer.time();
    return logp;
}


double SpimapModel::topologyPrior()
{
    Timer timer;
    double logp = birthDeathTreePriorFull(tree, stree, recon, events, 
                                          dupprob, lossprob,
                                          doomtable);
    top_runtime += timer.time();
    return logp;
}


//=============================================================================
// HKY sequence likelihood

HkySeqLikelihood::HkySeqLikelihood(int nseqs, int seqlen, char **seqs, 
                                   float *bgfreq, float tsvratio, int maxiter,
                                   double minlen, double maxlen) :
    nseqs(nseqs),
    seqlen(seqlen),
    seqs(seqs),
    bgfreq(bgfreq),
    tsvratio(tsvratio),
    maxiter(maxiter),
    minlen(minlen),
    maxlen(maxlen)
{}


double HkySeqLikelihood::findLengths(Tree *tree)
{ 
    return findMLBranchLengthsHky(tree, nseqs, seqs, bgfreq, 
                                  tsvratio, maxiter, minlen, maxlen);
}



} // namespace spidir

