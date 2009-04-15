#ifndef SPIDIR_BRANCH_PRIOR_H
#define SPIDIR_BRANCH_PRIOR_H

#include "Tree.h"
#include "spidir.h"
#include "birthdeath.h"

namespace spidir {

typedef void (*geneRateCallback) (float generate, Tree *tree, void *userdata);

extern "C" {



// fractional branches
enum {
    FRAC_NONE,
    FRAC_DIFF,
    FRAC_PARENT,
    FRAC_NODE
};



// Branch distribution parameters for one branch
class BranchParams
{
public:
    BranchParams(float _alpha=-1.0, float _beta=-1.0) :
        alpha(_alpha),
        beta(_beta)
    {}
    
    bool isNull()
    {
        return alpha == -1.0;
    }
    
    float alpha;
    float beta;
};


// Reconciliation parameters
class ReconParams
{
public:
    ReconParams(int nnodes, SpidirParams *params) :
        nnodes(nnodes),
	params(params),
        unfold(-1),
        unfolddist(0)
    {
        startspecies = new int [nnodes];
        midspecies = new ExtendArray<int> [nnodes];
        endspecies = new int [nnodes];
        
        startfrac = new int [nnodes];
        endfrac = new int [nnodes];
        midpoints = new float [nnodes];
        
        freebranches = new bool [nnodes];
    }
    
    ~ReconParams()
    {
        delete [] startspecies;
        delete [] midspecies;
        delete [] endspecies;
        
        delete [] startfrac;
        delete [] endfrac;
        delete [] midpoints;
        
        delete [] freebranches;
    }
    
    int nnodes;
    SpidirParams *params;

    int *startspecies;
    ExtendArray<int> *midspecies;
    int *endspecies;

    int *startfrac;
    int *endfrac;
    float *midpoints;

    bool *freebranches;
    int unfold;
    float unfolddist;
};




float branchPrior(int nnodes, int *ptree, float *dists,
		  int nsnodes, int *pstree, float *sdists,
		  int *recon, int *events,
		  float *sp_alpha, float *sp_beta, float generate, 
		  float predupprob=1.0, float dupprob=1.0, float lossprob=1.0,
		  float gene_alpha=0, float gene_beta=0);

} // extern "C"

float branchPrior(Tree *tree,
		  SpeciesTree *stree,
		  int *recon, int *events, SpidirParams *params,
		  float generate, 
		  float predupprob, float dupprob, float lossprob);

// get nodes in preorder (starting with given node)
//void getSubtree(int **ftree, int node, int *events, ExtendArray<int> *subnodes);

void getSubtree(Node *node, int *events, ExtendArray<Node*> *subnodes);

BranchParams getBranchParams(int node, Tree *tree, ReconParams *reconparams);

// Reconcile a branch to the species tree
void reconBranch(int node, 
		 Tree *tree, SpeciesTree *stree, 
		 int *recon, int *events, 
                 SpidirParams *params,
                 ReconParams *reconparams);

void determineFreeBranches(Tree *tree, SpeciesTree *stree, 
                           int *recon, int *events, float generate,
                           int *unfold, float *unfolddist, bool *freebranches);

void setRandomMidpoints(int root, Tree *tree,
                        Node **subnodes, int nsubnodes, 
                        int *recon, int *events, 
                        ReconParams *reconparams);


float rareEventsLikelihood(Tree *tree, SpeciesTree *stree, int *recon, 
                           int *events,
                           float predupprob, float dupprob, float lossprob);

float rareEventsLikelihood_old(Tree *tree, SpeciesTree *stree, int *recon, 
                               int *events,
                               float predupprob, float dupprob);

float maxPosteriorGeneRate(Tree *tree, SpeciesTree *stree,
                           int *recon, int *events, SpidirParams *params);

float maxPosteriorGeneRate(int nnodes, int *ptree, float *dists,
                           int nsnodes, int *pstree, 
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta);


void samplePosteriorGeneRate(Tree *tree,
                             int nseqs, char **seqs,
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *gene2species,
                             SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata);


void samplePosteriorGeneRate(Tree *tree,
                             int nseqs, char **seqs, 
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *recon, int *events, SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata);


void generateBranchLengths(int nnodes, int *ptree, 
                           int nsnodes, int *pstree,
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta,
                           float *dists);

void generateBranchLengths(Tree *tree,
                           SpeciesTree *stree,
                           int *recon, int *events,
                           SpidirParams *params,
                           float generate=-1.0, 
                           int subnode=-1, int subchild=-1);





} // namespace spidir

#endif // SPIDIR_BRANCH_PRIOR_H
