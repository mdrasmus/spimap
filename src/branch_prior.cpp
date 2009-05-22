//=============================================================================
//  SPIDIR - branch prior


//#include <gsl/gsl_integration.h>

// c++ headers
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>

// spidir headers
#include "common.h"
#include "branch_prior.h"
#include "birthdeath.h"
#include "phylogeny.h"
#include "ExtendArray.h"
#include "seq_likelihood.h"


namespace spidir {


BranchParams NULL_PARAM;




//=============================================================================
// gene rate estimation


// solves x^3 + ax^2 + bx + x = 0 for x
//   in our case there is only one positive root; find it
float maxCubicRoot(float a, float b, float c)
{
    float const accuracy = .001;
    float x = 0.0;
    float y;
    
    // first estimate should be negative    
    assert (x*x*x + a*x*x + b*x + c <= 0);
    
    // increase x until y is positive
    x = .01;
    do {
        x *= 2.0;
        y = x*x*x + a*x*x + b*x + c;
    } while (y < 0);
    
    // binary search to find root
    float min = x / 2.0;
    float max = x;
    
    while (max - min > accuracy) {
        x = (max + min) / 2.0;
        y = x*x*x + a*x*x + b*x + c;
        
        if (y == 0)
            return x;
        else if (y > 0)
            max = x;
        else
            min = x;
    }
    
    return x;
}

int floatcmp(const void *a, const void *b)
{
    float fa = *((float*)a);
    float fb = *((float*)b);
    
    if (fa < fb)
        return -1;
    else if (fa > fb)
        return 1;
    else
        return 0;
}


float mleGeneRate(int count, float *dists, float *means, float *sdevs,
                  float alpha, float beta)
{
    float a, b, c;
    float threshold = 0;
    
    a = (1.0 - alpha) / beta;

    // b = sum(means[i] * lens[i] / sdevs[i]**2
    //         for i in range(len(lens))) / beta    
    // c = - sum(lens[i] ** 2 / sdevs[i] ** 2
    //          for i in range(len(lens))) / beta    
    

    ExtendArray<float> dists2(0, count);
    dists2.extend(dists, count);
    qsort((void*) dists2.get(), dists2.size(), sizeof(float), floatcmp);
    //printFloatArray(dists2, dists2.size());
    //printFloatArray(dists, dists2.size());
    int limit = int(count * .5) + 1;
    if (limit < 4) limit = 4;
    threshold = dists2[limit];
    //printf("threshold %f\n", threshold);
    
    
    b = 0.0;
    c = 0.0;    
    for (int i=0; i<count; i++) {
        if (dists[i] > threshold && sdevs[i] > 0.0001) {
            b += means[i] * dists[i] / (sdevs[i] * sdevs[i]);
            c += dists[i] * dists[i] / (sdevs[i] * sdevs[i]);
        }
    }
    b /= beta;
    c = -c / beta;
    
    return maxCubicRoot(a, b, c);
}


// help set depths for every node
// depth is distance to next subtree root
void estimateGeneRate_helper(Tree *tree, Node *node, float *depths, int *sroots,
                             int *recon, int *events, bool free)
{
    int name = node->name;
    
    if (events[name] == EVENT_SPEC)
        free = false;
    
    if (node != tree->root) {
        int parent = node->parent->name;
        
        if (free) {
            // mark freebranches with -1
            depths[name] = -1;
            sroots[name] = sroots[parent];
        } else {
            if (events[parent] == EVENT_DUP) {
                depths[name] = depths[parent] + node->dist;
                sroots[name] = sroots[parent];
            } else {
                depths[name] = node->dist;
                sroots[name] = recon[parent];
            }
        }
    }
    
    for (int i=0; i<node->nchildren; i++)
        estimateGeneRate_helper(tree, node->children[i], 
                                depths, sroots, recon, events, free);    
}


float estimateGeneRate(Tree *tree, SpeciesTree *stree, 
                       int *recon, int *events, SpidirParams *params)
{
    float *depths = new float [tree->nnodes];
    int *sroots = new int [tree->nnodes];   // species roots
    
    depths[tree->root->name] = 0;
    sroots[tree->root->name] = recon[tree->root->name];
    
    // determine if top subtree is free
    bool free = (recon[tree->root->name] == stree->root->name && 
                 events[tree->root->name] == EVENT_DUP);
    
    estimateGeneRate_helper(tree, tree->root, depths, sroots, 
                            recon, events, free);
    
    //printFloatArray(depths, tree->nnodes);
    //printIntArray(sroots, tree->nnodes);
    
    float *dists = new float [tree->nnodes];
    float *means = new float [tree->nnodes];
    float *sdevs = new float [tree->nnodes];
    
    // make quick access to params
    float *mu = params->sp_alpha;
    float *sigma = params->sp_beta;
    
    int count = 0;
    
    for (int i=0; i<tree->nnodes; i++) {
        if (events[i] != EVENT_DUP && i != tree->root->name) {
            // we are at subtree leaf
            
            // figure out species branches that we cross
            // get total mean and variance of this path            
            float u = 0.0;
            float s2 = 0.0;   
            int snode = recon[i];
            
            // branch is also free if we do not cross any more species
            // don't estimate baserates from extra branches
            if (snode != sroots[i] && depths[i] != -1) {
                while (snode != sroots[i] && snode != stree->root->name) {
                    u += mu[snode];
                    s2 += sigma[snode]*sigma[snode];
                    snode = stree->nodes[snode]->parent->name;
                }
                assert(fabs(s2) > .0000001);
                                
                // save dist and params
                dists[count] = depths[i];
                means[count] = u;
                sdevs[count] = sqrt(s2);
                count++;
            }
        }
    }
    
    
    float generate = (count > 0) ? mleGeneRate(count, dists, means, sdevs, 
                                   params->gene_alpha, params->gene_beta) :
                                   0.0; // catch unusual case
    
    delete [] depths;
    delete [] sroots;
    delete [] dists;
    delete [] means;
    delete [] sdevs;    
    
    return generate;
}




//=============================================================================
// branch prior functions


float approxGammaSum(int nparams, float x, float *gs_alpha, float *gs_beta,
		     bool approx)
{
    const float tol = .001;
    const float minfrac = .01;

    // filter for extreme parameters
    float mean = 0.0;
    for (int i=0; i<nparams; i++) {
	if (isinf(gs_beta[i]) || isnan(gs_beta[i])) {
	    gs_alpha[i] = gs_alpha[--nparams];
	    gs_beta[i] = gs_beta[nparams];
	    i--;
	} else {
	    mean += gs_alpha[i] / gs_beta[i];
	}
    }

    // remove params with effectively zero mean
    float mu = 0.0;
    float var = 0.0;
    for (int i=0; i<nparams; i++) {	
	if (gs_alpha[i] / gs_beta[i] < minfrac * mean) {
	    gs_alpha[i] = gs_alpha[--nparams];
	    gs_beta[i] = gs_beta[nparams];
	    i--;
	} else {
	    mu += gs_alpha[i] / gs_beta[i];
	    var += gs_alpha[i] / gs_beta[i] / gs_beta[i];
	}
    }    
    

    // there is nothing to do
    if (nparams == 0)
	return -INFINITY;
    
    //float dist = node->dist / generate;
    
    float logp;
    if (approx) {
	// approximation
	float a2 = mean*mean/var;
	float b2 = mean/var;
	logp = gammalog(x, a2, b2);
    } else {
	logp = log(gammaSumPdf(x, nparams, gs_alpha, gs_beta, tol));
    }
    
    //float diff = fabs(logp - logp2);

    /*
    if (diff > log(1.3)) {
	printf("\n"
	       "n=%d; l=%f; g=%f\n", n, node->dist, generate);
	printf("t: "); printFloatArray(times, times.size());
	printf("a: "); printFloatArray(gs_alpha, nparams);
	printf("b: "); printFloatArray(gs_beta, nparams);
        
	printf("logp = %f\n", logp);
	printf("logp2 = %f\n", logp2);
	printf("diff = %f\n", diff);
	}*/

    if(isnan(logp))
	logp = -INFINITY;
    return logp;
}



// Generate a random sample of duplication points
void setRandomMidpoints(int root, Tree *tree, SpeciesTree *stree,
                        Node **subnodes, int nsubnodes, 
                        int *recon, int *events, 
                        ReconParams *reconparams,
			float birth, float death)
{
    const float esp = .0001;
    
    // deal with pre-duplications
    if (root == tree->root->name) {
	do {
	    reconparams->pretime = 
		expovariate(reconparams->params->pretime_lambda);
	} while (reconparams->pretime <= 0.0);
    }

    for (int i=0; i<nsubnodes; i++) {
	Node *node = subnodes[i];
        
        if (events[node->name] == EVENT_DUP) {
            float lastpoint;
            
            if (node->parent != NULL && 
		recon[node->name] == recon[node->parent->name])
                // if I'm the same species branch as my parent 
                // then he is my last midpoint
                lastpoint = reconparams->midpoints[node->parent->name];
            else
                // i'm the first on this branch so the last midpoint is zero
                lastpoint = 0.0;
            
            // pick a midpoint based on a DB process
            float remain = 1.0 - lastpoint;
	    int snode = recon[node->name];
	    float time;

	    if (snode == stree->root->name) {
		// use preduplication time for pre-duplicates
		time = reconparams->pretime;
	    } else {
		time = stree->nodes[snode]->dist;
	    }

	    assert(time > 0.0);

	    reconparams->midpoints[node->name] = lastpoint + esp * remain +
		sampleBirthWaitTime1(remain * time * (1.0 - esp), 
				     birth, death) / time;

        } else {
            // genes or speciations reconcile exactly to the end of the branch
            reconparams->midpoints[node->name] = 1.0;
        }
    }
}



// Reconcile a branch to the species tree
void reconBranch(int node, Tree *tree, SpeciesTree *stree,
		 int *recon, int *events, 
		 SpidirParams *params,
		 ReconParams *reconparams)
{
    Node** nodes = tree->nodes.get();
    Node** snodes = stree->nodes.get();
    
    ExtendArray<BranchPart> &parts = reconparams->parts[node];
    parts.clear();
    
    // set fractional branches
    if (recon[node] == recon[nodes[node]->parent->name]) {
        // start reconciles to a subportion of species branch
        if (events[node] == EVENT_DUP) {
            // only case k's are dependent
	    parts.append(BranchPart(recon[node], FRAC_DIFF));
        } else {
	    parts.append(BranchPart(recon[node], FRAC_PARENT));
	}

	// end reconcilies to nothing
    } else {
        if (events[nodes[node]->parent->name] == EVENT_DUP) {
            // start reconciles to last part of species branch
            int snode = recon[nodes[node]->parent->name];
	    parts.append(BranchPart(snode, FRAC_PARENT));
        }
	// else: start reconcilies to nothing

        if (events[node] == EVENT_DUP) {
            // end reconciles to first part of species branch
	    parts.append(BranchPart(recon[node], FRAC_NODE));
        }
	// else: end reconciles to nothing
    }
    
    // set midparams
    if (recon[node] != recon[nodes[node]->parent->name]) {
        // we begin and end on different branches
        int snode;

        // determine most recent species branch which we fully recon to
        if (events[node] == EVENT_DUP)
            snode = snodes[recon[node]]->parent->name;
        else
            snode = recon[node];

        // walk up species spath until starting species branch
        // starting species branch is either fractional or NULL
        int parent_snode;
        if (nodes[node]->parent != NULL)
            parent_snode = recon[nodes[node]->parent->name];
        else
            parent_snode = -1;
	
        while (snode != parent_snode && snodes[snode]->parent != NULL) {
	    parts.append(BranchPart(snode, FRAC_ONE));
            snode = snodes[snode]->parent->name;
        }
    } 

}



// get times vector from midpoints
void getReconTimes(Tree *tree, SpeciesTree *stree, 
		   Node *node, 
		   ReconParams *reconparams,
		   ExtendArray<float> &times)
{
    const float *k = reconparams->midpoints;
    const int name = node->name;
    const ExtendArray<BranchPart> &parts = reconparams->parts[name];


    times.clear();   

    for (int i=0; i<parts.size(); i++) {
	float time = stree->nodes[parts[i].species]->dist;
	
	if (parts[i].species == stree->root->name)
	    time = reconparams->pretime;

	switch (parts[i].frac) {
	case FRAC_DIFF:
	    times.append((k[name] - k[node->parent->name]) *
			 time);
	    break;

	case FRAC_PARENT:
	    float kp = k[node->parent->name];
	    times.append((1.0 - kp) * time);
	    break;

	case FRAC_NODE:
	    times.append(k[name] * time);
	    break;

	case FRAC_ONE:
	    times.append(time);
	    break;
	}
    }
}


// get gamma sum parameters
void getReconParams(Tree *tree, Node *node, 
		    ReconParams *reconparams,
		    float generate,
		    float *times,
		    float *gs_alpha,
		    float *gs_beta)
{
    SpidirParams *params = reconparams->params;
    ExtendArray<BranchPart> &parts = reconparams->parts[node->name];

    for (int j=0; j<parts.size(); j++) {
	int snode = parts[j].species;
	gs_alpha[j] = params->sp_alpha[snode];
	gs_beta[j] = (params->sp_beta[snode] / (generate * times[j]));
    }
}



// Calculate branch probability
float branchprob(Tree *tree, SpeciesTree *stree, Node *node,
		 float generate, ReconParams *reconparams,
		 bool approx=true)
{
    // get times
    ExtendArray<float> times(0, tree->nnodes);
    getReconTimes(tree, stree, node, reconparams, times);
    int nparams = times.size();

    // get gammaSum terms 
    float gs_alpha[nparams];
    float gs_beta[nparams];
    getReconParams(tree, node, reconparams,
		   generate, times, gs_alpha, gs_beta);
   
    // compute gamma sum
    return approxGammaSum(nparams, node->dist, gs_alpha, gs_beta, approx);
}


// Calculate branch probability
float branchprob_unfold(Tree *tree, SpeciesTree *stree, 
			float generate, ReconParams *reconparams,
			bool approx=true)
{
    Node *node0 = tree->root->children[0];
    Node *node1 = tree->root->children[1];

    // get times
    ExtendArray<float> times0(0, tree->nnodes);
    ExtendArray<float> times1(0, tree->nnodes);
    getReconTimes(tree, stree, node0, reconparams, times0);
    getReconTimes(tree, stree, node1, reconparams, times1);
    int nparams = times0.size() + times1.size();

    // get gamma sum terms 
    float gs_alpha[nparams];
    float gs_beta[nparams];
    getReconParams(tree, node0, reconparams,
		   generate, times0, gs_alpha, gs_beta);
    getReconParams(tree, node1, reconparams,
		   generate, times1, 
		   gs_alpha + times0.size(), 
		   gs_beta + times0.size());

    // compute gamma sum
    return approxGammaSum(nparams, node0->dist + node1->dist, 
			  gs_alpha, gs_beta, approx);
}


// get roots of speciation subtrees
void getSpecSubtrees(Tree *tree, int *events, ExtendArray<Node*> *rootnodes)
{
    for (int i=0; i<tree->nnodes; i++) {
	if (i == tree->root->name) {
	    rootnodes->append(tree->nodes[i]);
	} else if (events[i] == EVENT_SPEC) {
	    for (int j=0; j<2; j++) {
		rootnodes->append(tree->nodes[i]->children[j]);
	    }
	}
    }
}


// Returns True if duplications encounters
bool getSubtree(Node *node, int *events, ExtendArray<Node*> *subnodes)
{
    subnodes->append(node);
    bool dupsPresent = false;

    // recurse
    if (events[node->name] == EVENT_DUP || node->parent == NULL) {
	if (events[node->name] == EVENT_DUP)
	    dupsPresent = true;
	
        for (int i=0; i<node->nchildren; i++) 
            getSubtree(node->children[i], events, subnodes);
    }

    return dupsPresent;
}


class BranchPriorCalculator
{
public:
    BranchPriorCalculator(Tree *_tree,
			  SpeciesTree *_stree,
			  int *_recon, int *_events, 
			  SpidirParams *_params,
			  float _birth, float _death,
			  int _nsamples=1000, bool _approx=true) :
        tree(_tree),
        stree(_stree),
        recon(_recon),
        events(_events),
        params(_params),
	birth(_birth),
	death(_death),
	reconparams(_tree->nnodes, _params),
	rootnodes(0, _tree->nnodes),
	nsamples(_nsamples),
	approx(_approx)
    {
	// determine speciation subtrees	
	getSpecSubtrees(tree, events, &rootnodes);
    }
    
    
    ~BranchPriorCalculator()
    {
    }
    


    // subtree prior conditioned on divergence times
    float subtreeprior_cond(Tree *tree, SpeciesTree *stree,
			    int *recon, float generate,
			    ReconParams *reconparams,
			    ExtendArray<Node*> &subnodes,
			    bool unfold)
    {
	float logp = 0.0;
        
	// loop through all branches in subtree
	for (int j=0; j<subnodes.size(); j++) {
	    Node *node = subnodes[j];
	
	    //if (recon[node->name] == sroot)
	    //	continue;

	    if (node == tree->root)
		continue;

	    if (unfold) {
		if (node == tree->root->children[1]) {
		    // skip right branch
		    continue;
		} else if (node == tree->root->children[0]) {
		    // unfold left branch
		    logp += branchprob_unfold(tree, stree,
					      generate, reconparams);
		    continue;
		}
	    }
	    logp += branchprob(tree, stree, node, generate, reconparams);
	}
            

	return logp;
    }


    // Calculate the likelihood of a subtree
    float subtreeprior(int root, float generate, ReconParams *reconparams)
    {
   
	// set reconparams by traversing subtree
	ExtendArray<Node*> subnodes(0, tree->nnodes);
	bool dupsPresent = getSubtree(tree->nodes[root], events, &subnodes);
    
	// reconcile each branch
	for (int i=0; i<subnodes.size(); i++)
	    if (subnodes[i] != tree->root)
		reconBranch(subnodes[i]->name, tree, stree, 
			    recon, events, params, reconparams);
	
	// root branch must be unfolded
	bool unfold = (root == tree->root->name &&
		       tree->root->nchildren == 2);

	
	if (!dupsPresent) {
	    return subtreeprior_cond(tree, stree, recon, 
				     generate, reconparams, subnodes,
				     unfold);
	} else {
        
	    // choose number of samples based on number of nodes 
	    // to integrate over
	    //nsamples = int(500*logf(subnodes.size())) + 200;
	    //if (nsamples > 2000) nsamples = 2000;
        	
	    // perform integration by sampling
	    double prob = 0.0;
	    for (int i=0; i<nsamples; i++) {
		double sampleLogl = 0.0;
            
		// propose a setting of midpoints
		setRandomMidpoints(root, tree, stree,
				   subnodes, subnodes.size(),
				   recon, events, reconparams,
				   birth, death);
            

		sampleLogl = subtreeprior_cond(tree, stree, recon, 
					       generate, reconparams, 
					       subnodes,
					       unfold);
	    
		prob += exp(sampleLogl) / nsamples;
	    }
        
	    return log(prob);
	}
    }


    
    double calc_cond(float generate)
    {
        double logp = 0.0;
	
	// multiple the probability of each subtree
	for (int i=0; i<rootnodes.size(); i++) {
	    logp += subtreeprior(rootnodes[i]->name, generate, &reconparams);
	}
        
        // gene rate probability
        if (params->gene_alpha > 0 && params->gene_beta > 0)
            logp += gammalog(generate, params->gene_alpha, params->gene_beta);
	
        printLog(LOG_HIGH, "generate: %f %f\n", generate, exp(logp));
        return logp;
    }
    
    inline double calcprob_cond(float generate)
    {
        return exp(calc_cond(generate));
    }


    double calc()
    {
	double maxprob = -INFINITY;
	float argmax_generate = params->gene_alpha / params->gene_beta;
	float est_generate = estimateGeneRate(tree, stree, recon, events, params);
        printLog(LOG_HIGH, "est_generate: %f\n", est_generate);
        
	// TODO: make this equal portions of gamma
        double logp = -INFINITY;       
        float gstart = params->gene_alpha / params->gene_beta * 0.05;
        float gend = params->gene_alpha / params->gene_beta * 3.0;
        float step = (gend - gstart) / 20.0;


	// integrate over gene rates
        for (float g=gstart; g<gend; g+=step) {
            float gi = g + step / 2.0;
            double p = calc_cond(gi);
	    
            logp = logadd(logp, p);
            printLog(LOG_HIGH, "generate_int: %f %f\n", gi, p);
            if (p > maxprob) {
                maxprob = p;
                argmax_generate = gi;
            }
        }
        
        // multiply probabilty by integration step
        logp += log(step);
        
        printLog(LOG_HIGH, "argmax gene rate: %f\n", argmax_generate);
        if (est_generate > 0.0)
            printLog(LOG_HIGH, "branch prior: %f %f %f\n", maxprob, logp,
		     calc_cond(est_generate));

	return logp;
    }
    
    
protected:
    Tree *tree;    
    SpeciesTree *stree;
    int *recon;
    int *events;
    SpidirParams *params;  
    float birth;
    float death;    
    ReconParams reconparams;
    ExtendArray<Node*> rootnodes;

    int nsamples;
    bool approx;
  
};

  

// TODO: move birth, death into param
float branchPrior(Tree *tree,
		  SpeciesTree *stree,
		  int *recon, int *events, SpidirParams *params,
		  float generate,
		  float pretime_lambda, float birth, float death,
		  int nsamples, bool approx)
{
    float logl = 0.0; // log likelihood
    Timer timer;

    
    BranchPriorCalculator priorcalc(tree, stree, recon, events, params,
				    birth, death, nsamples, approx);
    
    if (generate > 0) {
        // estimate with given generate
        logl = priorcalc.calc_cond(generate);
    } else {
        logl = priorcalc.calc();
    }
    
    //setLogLevel(LOG_MEDIUM);
    printLog(LOG_MEDIUM, "branch prior time: %f\n", timer.time());
    
    return logl;
}


extern "C" {

// Calculate the likelihood of a tree
float branchPrior(int nnodes, int *ptree, float *dists,
		  int nsnodes, int *pstree, float *sdists,
		  int *recon, int *events,
		  float *sp_alpha, float *sp_beta, float generate,
		  float pretime_lambda, float birth, float death,
		  float gene_alpha, float gene_beta,
		  int nsamples, bool approx)
{
    // create tree objects
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);
    
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();
    stree.setDists(sdists);
    
    SpidirParams params(nsnodes, NULL, sp_alpha, sp_beta, 
			gene_alpha, gene_beta, pretime_lambda);
    

    return branchPrior(&tree, &stree,
		       recon, events, &params, 
		       generate, 
		       pretime_lambda, birth, death,
		       nsamples, approx);
}

} // extern "C"









//=============================================================================
// Gene rate estimation


// functor of gene rate posterior derivative
class GeneRateDerivative 
{
public:
    GeneRateDerivative(Tree *tree, SpeciesTree *stree,
                       int *recon, int *events, SpidirParams *params) :
        lkcalc(tree, stree, recon, events, params, 0.0001, 0.0002)
    {
    }
    
    float operator()(float generate)
    {
        const float step = .05;
        return lkcalc.calc_cond(generate + step) - lkcalc.calc_cond(generate);
    }
    
    BranchPriorCalculator lkcalc;
};


// Ignores sequence data, assumes given branch lengths are perfectly known
float maxPosteriorGeneRate(Tree *tree, SpeciesTree *stree,
                           int *recon, int *events, SpidirParams *params)
{
    // get initial estimate of gene rate and a bound to search around
    float est_generate = estimateGeneRate(tree, stree, recon, events, params);
    float maxg = est_generate * 1.5; //params->alpha / params->beta * 2.0;
    float ming = est_generate / 1.5;
    
    // find actual max posterior of gene rate
    GeneRateDerivative df(tree, stree, recon, events, params);
    return bisectRoot(df, ming, maxg, (maxg - ming) / 1000.0);
}


float maxPosteriorGeneRate(int nnodes, int *ptree, float *dists,
                           int nsnodes, int *pstree, 
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta)
{
    // create tree objects
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);
    
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();    

    SpidirParams params = SpidirParams(nsnodes, NULL, mu, sigma, alpha, beta);
    
    return maxPosteriorGeneRate(&tree, &stree, recon, events, &params);
}


void samplePosteriorGeneRate(Tree *tree, 
                             int nseqs, char **seqs, 
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree, 
                             int *gene2species,
                             SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata)
{
    // use parsimonious reconciliation as default
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);

    // check events
    for (int i=0; i<tree->nnodes; i++) {
	assert(events[i] <= 2 && events[i] >= 0);
    }

    samplePosteriorGeneRate(tree,
                            nseqs, seqs, 
                            bgfreq, ratio,
                            stree,
                            recon, events, params,
                            nsamples,
                            callback,
                            userdata);
}


// Uses MCMC to sample from P(B,G|T,D)
void samplePosteriorGeneRate(Tree *tree,
                             int nseqs, char **seqs, 
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *recon, int *events, SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata)
{
    // state variables B, G
    ExtendArray<float> dists(tree->nnodes);
    float next_generate = 0, generate = 0;
    float logl = -INFINITY;
    float logl2 = 999;
    float next_logl, next_logl2;
    float alpha;

    // TODO: revive
    assert(0);

    const float generate_step = .2;
    const float min_generate = .0001;

    // initialize first state
    generate = gammavariate(params->gene_alpha, params->gene_beta);
    // TODO: replace
    //generateBranchLengths(tree, stree,
    //                      recon, events,
    //                      params, generate);

    
    // perform MCMC
    for (int i=0; i<nsamples; i++) {
        // generate a new state 
	//printf("sample %d\n", i);
     
	for (int j=0; j<tree->nnodes; j++) {
	    assert(events[j] <= 2 && events[j] >= 0);

	    // NAN distances
	    if (isnan(tree->nodes[j]->dist)) {
		tree->nodes[j]->dist = min_generate;
	    }	    
	}
 
        if (frand() < .5) {
	    //printf("sample G\n");
            printLog(LOG_HIGH, "sample G: ");

            // sample G_2
            next_generate = frand(max(generate - generate_step, min_generate),
                                  generate + generate_step);
            //next_generate = gammavariate(params->alpha, params->beta);

            // if P(B_1|T,G_1) not exist, make it
            if (logl2 > 0) {                
                logl2 = branchPrior(tree, stree, recon, events, params,
				    generate,
				    1, 1, 1);
            }
	    
            // set new branch lengths B_2
            for (int j=0; j<tree->nnodes; j++) {
                tree->nodes[j]->dist *= next_generate / generate;

		// sometimes when generate change is drastic we need to catch
		// NAN distances
		if (isnan(tree->nodes[j]->dist)) {
		    tree->nodes[j]->dist = min_generate;
		}
	    }

            
            // calculate P(D|B_2,T)
            next_logl = calcSeqProbHky(tree, nseqs, seqs, bgfreq, ratio);
	    
            // calculate P(B_2|T,G_2)
            next_logl2 = branchPrior(tree, stree, recon, events, params,
				     next_generate,
				     1, 1, 1);

            // q(G_1|G_2) = 1 /(next_generate + generate_step - 
            //                  max(next_generate - generate_step, 
            //                      min_generate))
            // q(G_2|G_1) = 1 /(generate + generate_step - 
            //                  max(generate - generate_step, 
            //                      min_generate))
            // qratio = (generate + generate_step - 
            //                  max(generate - generate_step, 
            //                      min_generate)) /
            //          (next_generate + generate_step - 
            //                  max(next_generate - generate_step, 
            //                      min_generate))


            float logl3 = gammalog(next_generate, 
                                   params->gene_alpha, params->gene_beta) - 
                          gammalog(generate, 
                                   params->gene_alpha, params->gene_beta) +
                          log((generate + generate_step - 
                               max(generate - generate_step, 
                                   min_generate)) /
                              (next_generate + generate_step - 
                               max(next_generate - generate_step, 
                                   min_generate)));

            alpha = exp(next_logl + next_logl2 - logl - logl2 + logl3);

        } else {
	    //printf("sample B\n");
            printLog(LOG_HIGH, "sample B: ");

            // keep same gene rate G_2 = G_1
            next_generate = generate;
            
            // sample B_2
	    // TODO: replace
            //generateBranchLengths(tree, stree,
            //                      recon, events,
            //                      params, next_generate, -2, -2);
	    
	    for (int j=0; j<tree->nnodes; j++) {
		// NAN distances
		if (isnan(tree->nodes[j]->dist)) {
		    tree->nodes[j]->dist = min_generate;
		}	    
	    }

            // calculate P(D|B_2,T)
            next_logl = calcSeqProbHky(tree, nseqs, seqs, bgfreq, ratio);
            next_logl2 = 999;
            
            alpha = exp(next_logl - logl);
        }
        

        if (frand() <= alpha) {
            // accept: logl, G, B
            printLog(LOG_HIGH, "accept\n");
            logl = next_logl;
            logl2 = next_logl2;
            generate = next_generate;
            tree->getDists(dists);
        } else {
            // reject
            printLog(LOG_HIGH, "reject\n");
            
            // restore previous B
            tree->setDists(dists);
        }
        
        // return sample
        callback(generate, tree, userdata);
    }
}




// Uses MCMC to sample from P(B,G|T,D)
void samplePosteriorGeneRate_old(Tree *tree,
                             int nseqs, char **seqs, 
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *recon, int *events, SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata)
{
    // state variables B, G
    ExtendArray<float> dists(tree->nnodes);
    float next_generate = 0, generate = 0;
    float logl = -INFINITY;
    float next_logl;
    
    // TODO: revive
    assert(0);

    for (int i=0; i<nsamples; i++) {
        // generate a new state 

        // sample G
        next_generate = gammavariate(params->gene_alpha, params->gene_beta);

        // sample B
	// TODO: replace
        //generateBranchLengths(tree, stree,
        //                      recon, events,
        //                      params, next_generate);
        
        // calculate P(D|B,T)
        next_logl = calcSeqProbHky(tree, nseqs, seqs, bgfreq, ratio);

        //printf(">> %f %f\n", next_logl, logl);

        float alpha = exp(next_logl - logl);
        if (frand() <= alpha) {
            // accept: logl, G, B
            printLog(LOG_HIGH, "accept: %f, %f\n", next_logl, logl);
            logl = next_logl;
            generate = next_generate;
            tree->getDists(dists);
        } else {
            // reject
            printLog(LOG_HIGH, "reject: %f, %f\n", next_logl, logl);
            
            // restore previous B
            tree->setDists(dists);
        }
        
        // return sample
        callback(generate, tree, userdata);
    }
}





} // namespace spidir
