
#ifndef SPIDIR_SEARCH_H
#define SPIDIR_SEARCH_H

#include "spidir.h"
#include <set>


namespace spidir {

using namespace std;


struct lttree
{
    bool operator()(const int* k1, const int* k2) const
    {
        for (int i=0; k1[i] != -2; i++) {
            if (k1[i] < k2[i])
                return true;
            if (k1[i] > k2[i])
                return false;
        }
        return false;
    }
};

class TreeSet
{
public:
    TreeSet() : key(100) {}
    ~TreeSet();

    void clear();
    bool insert(Tree *tree);
    bool has(Tree *tree);
    int size() { return trees.size(); }

    ExtendArray<int> key;
    typedef set<int*, lttree> Set;
    Set trees;
};


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
    virtual void testCorrect(Tree *tree)
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
    NniProposer(int niter=500);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more() { return iter < niter; }
    virtual void reset() { iter = 0; }


protected:    
    int niter;
    int iter;

    Node *nodea;
    Node *nodeb;
};


class SprProposer: public NniProposer
{
public:
    SprProposer(int niter=500);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);

protected:    
    int niter;
    int iter;

    Node *nodea;
    Node *nodeb;
    Node *nodec;
};


class SprNbrProposer: public NniProposer
{
public:
    SprNbrProposer(int niter=500, int radius=4);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    void reset() { 
        iter = 0; 
        reverted = false;
        basetree = NULL;
    }

    void pickNewSubtree();

protected:
    int radius;
    Tree *basetree;
    Node *subtree;
    list<Node*> queue;
    vector<int> pathdists;
    bool reverted;
};


class MixProposer: public TopologyProposer
{
public:
    MixProposer(int niter=500) : totalWeight(0), niter(niter), iter(0) {}

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);    
    virtual bool more() { return iter < niter; }
    virtual void reset() { 
        iter = 0; 

        // propagate reset
        for (unsigned int i=0; i<methods.size(); i++)
            methods[i].first->reset();
    }

    void addProposer(TopologyProposer *proposer, float weight);

protected:

    float totalWeight;
    typedef pair<TopologyProposer*,float> Method;
    vector<Method> methods;
    int lastPropose;
    int niter;
    int iter;
};


class ReconRootProposer: public TopologyProposer
{
public:
    ReconRootProposer(TopologyProposer *proposer,
                      SpeciesTree *stree=NULL, int *gene2species=NULL) :
        proposer(proposer),
        stree(stree),
        gene2species(gene2species),
        oldroot1(NULL),
        oldroot2(NULL)
    {}

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more() { return proposer->more(); }
    virtual void reset() { return proposer->reset(); }
    virtual void accept(bool accepted) { proposer->accept(accepted); }

protected:
    TopologyProposer *proposer;
    SpeciesTree *stree;
    int *gene2species;
    Node *oldroot1;
    Node *oldroot2;
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

    virtual ~DupLossProposer();

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more() { return iter < niter; }
    virtual void reset() {
        iter = 0;
        uniques.clear();
        proposer->reset();
    }
    
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
    double *doomtable;
    const static int maxdoom = 10;
    TreeSet uniques;

    ExtendArray<int> recon;
    ExtendArray<int> events;
    Tree *oldtop;
};


class UniqueProposer: public TopologyProposer
{
public:
    UniqueProposer(TopologyProposer *proposer, int niter=-1, int ntries=10) : 
        proposer(proposer), niter(niter), iter(0), ntries(ntries) {}
    virtual ~UniqueProposer();

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree) { return proposer->revert(tree); }
    virtual bool more() { 
        if (niter == -1)
            return proposer->more(); 
        else 
            return iter < niter;
    }
    virtual void reset() { 
        seenTrees.clear();
        proposer->reset(); 
    }
    
    TopologyProposer *proposer;
    TreeSet seenTrees;
    int niter;
    int iter;
    int ntries;
};

/*=============================================================================
NEW DUP/LOSS PROPOSER

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


=============================================================================*/


class BranchLengthFitter
{
public:
    BranchLengthFitter() : runtime(0.0) {}
    virtual ~BranchLengthFitter() {}
    virtual double findLengths(Tree *tree) {return 0.0;}

    float runtime;
};


class HkyFitter : public BranchLengthFitter
{
public:
    HkyFitter(int nseqs, int seqlen, char **seqs, 
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



//=============================================================================



class Prior
{
public:
    Prior() :
        branch_runtime(0),
        top_runtime(0)
    {}
    virtual ~Prior() {}
    
    virtual double branchPrior(Tree *tree) { return 0.0; }
    virtual double topologyPrior(Tree *tree) { return 0.0; }

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
		int nsample, bool approx, bool useBranchPrior);

    virtual ~SpidirPrior();

    virtual double branchPrior(Tree *tree);
    virtual double topologyPrior(Tree *tree);

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
    bool useBranchPrior;
    double *doomtable;
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

    TreeSearch() :
        proposal_runtime(0.0)
    {}

    virtual ~TreeSearch()
    {}

    virtual Tree *search(Tree *initTree, 
			 string *genes, 
			 int nseqs, int seqlen, char **seqs)
    { return NULL; }

    double proposal_runtime;
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
