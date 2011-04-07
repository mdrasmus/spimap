/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree evolution model

=============================================================================*/


#ifndef SPIDIR_MODEL_H
#define SPIDIR_MODEL_H

#include "model_params.h"
#include "newick.h"
#include <set>


namespace spidir {

using namespace std;



class SeqLikelihood
{
public:
    SeqLikelihood() {}
    virtual ~SeqLikelihood() {}
    virtual double findLengths(Tree *tree) {return 0.0;}
};


class HkySeqLikelihood : public SeqLikelihood
{
public:
    HkySeqLikelihood(int nseqs, int seqlen, char **seqs, 
                     float *bgfreq, float tsvratio, int maxiter, 
                     double minlen=0.0001, double maxlen=10.0);
    virtual double findLengths(Tree *tree);

    int nseqs;
    int seqlen;
    char **seqs;    
    float *bgfreq;
    float tsvratio;
    int maxiter;
    double minlen;
    double maxlen;
};




class Model
{
public:
    Model() :
        seq_runtime(0),
        branch_runtime(0),
        top_runtime(0)
    {}
    virtual ~Model() {}
    
    virtual void setTree(Tree *_tree) { tree = _tree; }
    virtual double likelihood() { return 0.0; }
    virtual double branchPrior() { return 0.0; }
    virtual double topologyPrior() { return 0.0; }

    virtual SpeciesTree *getSpeciesTree() { return NULL; }
    virtual int *getGene2species() { return NULL; }

    Tree *tree;

    // runtimes
    float seq_runtime;
    float branch_runtime;
    float top_runtime;
};


class SpimapModel : public Model
{
public:
    SpimapModel(int nnodes, SpeciesTree *stree, 
		SpidirParams *params, 
		int *gene2species,
		float predupprob, float dupprob, float lossprob,
		int nsample, bool approx, bool useBranchPrior);

    virtual ~SpimapModel();

    virtual void setTree(Tree *_tree);
    virtual double likelihood();
    virtual double branchPrior();
    virtual double topologyPrior();
    
    virtual SpeciesTree *getSpeciesTree() { return stree; }
    virtual int *getGene2species() { return gene2species; }

    void setLikelihoodFunc(SeqLikelihood *l) { likelihoodFunc = l; }
    
protected:
    int nnodes;
    SpeciesTree *stree;
    SpidirParams *params;
    int *gene2species;
    ExtendArray<int> recon;
    ExtendArray<int> events;
    float predupprob;
    float dupprob;
    float lossprob;
    int nsamples;
    bool approx;
    bool useBranchPrior;
    double *doomtable;
    SeqLikelihood *likelihoodFunc;
};


} // namespace

#endif // SPIDIR_MODEL_H
