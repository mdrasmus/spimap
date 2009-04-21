#ifndef SPIDIR_SEARCH_H
#define SPIDIR_SEARCH_H

#include "spidir.h"


namespace spidir {

class TopologyProposer
{
public:
    TopologyProposer() {}
    virtual ~TopologyProposer() {}
    virtual void propose(Tree *tree) {}
    virtual void revert(Tree *tree) {}
    virtual bool more() { return false; }
    virtual void reset() {}
    virtual void setCorrect(Tree *tree) {}
    virtual bool seenCorrect() { return false; }
};


class NniProposer: public TopologyProposer
{
public:
    NniProposer(SpeciesTree *stree=NULL, int *gene2species=NULL, int niter=500);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more();
    virtual void setCorrect(Tree *tree) { correctTree = tree; }
    virtual bool seenCorrect() { return correctSeen; }
    virtual void reset() { iter = 0; }

protected:    
    Node *nodea;
    Node *nodeb;
    Node *nodec;
    Node *noded;
    Node *oldroot1;
    Node *oldroot2;
    SpeciesTree *stree;
    int *gene2species;
    int niter;
    int iter;
    Tree *correctTree;
    bool correctSeen;
};


class SprNniProposer: public NniProposer
{
public:
    SprNniProposer(SpeciesTree *stree=NULL, int *gene2species=NULL, 
                   int niter=500, float sprRatio=0.5);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);    
    
protected:
    typedef enum {
        PROPOSE_NONE,
        PROPOSE_NNI,
        PROPOSE_SPR
    } ProposeType;
    
    float sprRatio;
    ProposeType lastPropose;
};


class DupLossProposer: public TopologyProposer
{
public:
    DupLossProposer(TopologyProposer *proposer, 
                    SpeciesTree *stree,
                    int *gene2species,
                    float dupprob,
                    float lossprob,
                    int quickiter=100, int niter=500);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual void setCorrect(Tree *tree) { correctTree = tree; }
    virtual bool seenCorrect() { return correctSeen; }    
    virtual bool more() { return iter < niter; }
    virtual void reset() { iter = 0; }    
    
protected:
    TopologyProposer *proposer;
    int quickiter;
    int niter;
    int iter;
    Tree *correctTree;
    bool correctSeen;
    SpeciesTree *stree;
    int *gene2species;
    float dupprob;
    float lossprob;

    ExtendArray<int> recon;
    ExtendArray<int> events;
    Tree *oldtop;
};


class BranchLengthFitter
{
public:
    BranchLengthFitter() {}
    virtual ~BranchLengthFitter() {}
    virtual float findLengths(Tree *tree) {return 0.0;}
};


class ParsimonyFitter : public BranchLengthFitter
{
public:
    ParsimonyFitter(int nseqs, int seqlen, char **seqs);
    virtual float findLengths(Tree *tree);
    
    int nseqs;
    int seqlen;
    char **seqs;
};


class HkyFitter : public BranchLengthFitter
{
public:
    HkyFitter(int nseqs, int seqlen, char **seqs, 
              float *bgfreq, float tsvratio, int maxiter, bool useLogl=true);
    virtual float findLengths(Tree *tree);

    int nseqs;
    int seqlen;
    char **seqs;    
    float *bgfreq;
    float tsvratio;
    int maxiter;
    bool useLogl;
};


class SpidirSample : public BranchLengthFitter
{
public:
    SpidirSample(SpeciesTree *stree, SpidirParams *params, int *gene2species) :
        stree(stree),
        params(params),
        gene2species(gene2species)
    {}
    virtual float findLengths(Tree *tree);
    
    SpeciesTree *stree;
    SpidirParams *params;
    int *gene2species;
};


class HkySpidirSample : public BranchLengthFitter
{
public:
    HkySpidirSample(SpeciesTree *stree, SpidirParams *params, int *gene2species,
                    int nseqs, int seqlen, char **seqs, 
                    float *bgfreq, float tsvratio, int maxiter) :
        stree(stree),
        params(params),
        gene2species(gene2species),
        nseqs(nseqs),
        seqlen(seqlen),
        seqs(seqs),
        bgfreq(bgfreq),
        maxiter(maxiter)
        
    {}
    virtual float findLengths(Tree *tree);
    
    SpeciesTree *stree;
    SpidirParams *params;
    int *gene2species;
    int nseqs;
    int seqlen;
    char **seqs;    
    float *bgfreq;
    float tsvratio;
    int maxiter;
};


class BirthDeathFitter : public BranchLengthFitter
{
public:
    BirthDeathFitter(int nseqs, int seqlen, char **seqs, 
                     float *bgfreq, float tsvratio,
                     SpeciesTree *stree, int *gene2species,
                     float birthRate, float deathRate);
    virtual float findLengths(Tree *tree);
    
    int nseqs;
    int seqlen;
    char **seqs;    
    float *bgfreq;
    float tsvratio;    
    SpeciesTree *stree;
    int *gene2species;
    float birthRate;
    float deathRate;
};


//=============================================================================


class Prior
{
public:
    Prior() {}
    virtual ~Prior() {}
    
    virtual float branchPrior(Tree *tree) { return 0.0; }
    virtual float topologyPrior(Tree *tree) { return 0.0; }

    virtual SpeciesTree *getSpeciesTree() { return NULL; }
    virtual int *getGene2species() { return NULL; }    
};


class SpidirPrior : public Prior
{
public:
    SpidirPrior(int nnodes, SpeciesTree *stree, 
		SpidirParams *params, 
		int *gene2species,
		float predupprob, float dupprob, float lossprob,
		bool estGenerate);

    virtual ~SpidirPrior();

    virtual float branchPrior(Tree *tree);
    virtual float topologyPrior(Tree *tree);

    virtual SpeciesTree *getSpeciesTree() { return stree; }
    virtual int *getGene2species() { return gene2species; }
    
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
    bool estGenerate;
    float *doomtable;
};


/*
class HkyBranchLikelihoodFunc : public BranchLikelihoodFunc
{
public:
    HkyBranchLikelihoodFunc(int nseqs, int seqlen, char **seqs, 
                         float *bgfreq, float tsvratio) :
        nseqs(nseqs),
        seqlen(seqlen),
        seqs(seqs),
        bgfreq(bgfreq),
        tsvratio(tsvratio)
    {}
    
    virtual float likelihood(Tree *tree);
    virtual float likelihood2(Tree *tree) { return 0.0; }
    virtual SpeciesTree *getSpeciesTree() { return NULL; }
    virtual int *getGene2species() { return NULL; }

    int nseqs;
    int seqlen;
    char **seqs;    
    float *bgfreq;
    float tsvratio; 
};

*/

class SampleFunc
{
public:
    SampleFunc(FILE *output) :
        output(output)
    {
    }
    
    virtual ~SampleFunc()
    {
        fclose(output);
    }

    void operator()(Tree *tree)
    {
        tree->writeNewick(output, NULL, 0, true);
        fprintf(output, "\n");
    }
    
protected:
    FILE *output;
};


class TreeSearch
{
public:

    TreeSearch()
    {}

    virtual ~TreeSearch()
    {}

    virtual Tree *search(Tree *initTree, 
			 string *genes, 
			 int nseqs, int seqlen, char **seqs)
    { return NULL; }

};


class TreeSearchClimb : public TreeSearch
{
public:

    TreeSearchClimb(Prior *prior,
		    TopologyProposer *proposer,
		    BranchLengthFitter *fitter);

    virtual ~TreeSearchClimb();

    virtual Tree *search(Tree *initTree, 
			 string *genes, 
			 int nseqs, int seqlen, char **seqs);


protected:
    Prior *prior;
    TopologyProposer *proposer;
    BranchLengthFitter *fitter;

};




Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs,
                     SpeciesTree *stree, int *gene2species);
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs);




Tree *searchMCMC(Tree *initTree, 
                 string *genes, int nseqs, int seqlen, char **seqs,
                 SampleFunc *samples,
                 Prior *lkfunc,
                 TopologyProposer *proposer,
                 BranchLengthFitter *fitter);


} // namespace spidir

#endif // SPIDIR_SEARCH_H
