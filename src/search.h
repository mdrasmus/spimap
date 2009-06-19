#ifndef SPIDIR_SEARCH_H
#define SPIDIR_SEARCH_H

#include "spidir.h"
#include <list>
#include <vector>


namespace spidir {

using namespace std;


//=============================================================================
// proporsers

class TopologyProposer
{
public:
    TopologyProposer() :
        correctTree(NULL),
        correctSeen(false)
    {}

    virtual ~TopologyProposer() {}
    virtual void propose(Tree *tree) {}
    virtual void revert(Tree *tree) {}
    virtual bool more() { return false; }
    virtual void reset() {}
    virtual void accept(bool accepted) {}

    virtual void setCorrect(Tree *tree) { correctTree = tree; }
    virtual Tree *getCorrect() { return correctTree; }
    virtual bool seenCorrect() { return correctSeen; }
    virtual inline void testCorrect(Tree *tree)
    {
        // debug: keep track of correct tree in search
        if (correctTree) {
            if (tree->sameTopology(correctTree))
                correctSeen = true;
        }
    }

protected:
    Tree *correctTree;
    bool correctSeen;
};


class NniProposer: public TopologyProposer
{
public:
    NniProposer(SpeciesTree *stree=NULL, int *gene2species=NULL, int niter=500);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more();
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
};


class SprProposer: public NniProposer
{
public:
    SprProposer(SpeciesTree *stree=NULL, int *gene2species=NULL, 
                int niter=500);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
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



class SprNbrProposer: public NniProposer
{
public:
    SprNbrProposer(SpeciesTree *stree=NULL, int *gene2species=NULL, 
                   int niter=500, int radius=4);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual void reset();

    void pickNewSubtree();

protected:
    int radius;
    Tree *basetree;
    Node *subtree;
    list<Node*> queue;
    vector<int> pathdists;
};


class DupLossProposer: public TopologyProposer
{
public:
    DupLossProposer(TopologyProposer *proposer, 
                    SpeciesTree *stree,
                    int *gene2species,
                    float dupprob,
                    float lossprob,
                    int quickiter=100, int niter=500, int nsamples=10,
                    bool extend=true);
    virtual ~DupLossProposer();


    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more() { return iter < niter; }
    virtual void reset();
    virtual void accept(bool accepted);

    void queueTrees(Tree *tree);
    void clearQueue();

    typedef pair<Tree*,float> TreeProp;


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
    float *doomtable;
    const static int maxdoom = 10;
    vector<TreeProp> queue;
    float sum;
    ExtendArray<int> recon;
    ExtendArray<int> events;
    Tree *oldtop;
    int nsamples;
    int samplei;
    int treesize;
    bool extend;
};


//=============================================================================
// fitters

class BranchLengthFitter
{
public:
    BranchLengthFitter() :
        runtime(0)
    {}
    virtual ~BranchLengthFitter() {}
    virtual float findLengths(Tree *tree) { return 0.0; }

    float runtime;
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



//=============================================================================
// priors

class Prior
{
public:
    Prior() :
        branch_runtime(0),
        top_runtime(0)
    {}
    virtual ~Prior() {}
    
    virtual float branchPrior(Tree *tree) { return 0.0; }
    virtual float topologyPrior(Tree *tree) { return 0.0; }

    virtual SpeciesTree *getSpeciesTree() { return NULL; }
    virtual int *getGene2species() { return NULL; }


    float branch_runtime;
    float top_runtime;
};


class SpidirPrior : public Prior
{
public:
    SpidirPrior(int nnodes, SpeciesTree *stree, 
		SpidirParams *params, 
		int *gene2species,
		float predupprob, float dupprob, float lossprob,
		int nsample, bool approx, bool estGenerate);

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
    int nsamples;
    bool approx;
    bool estGenerate;
    float *doomtable;
};



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
