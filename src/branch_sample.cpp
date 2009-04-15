//=============================================================================
//  SPIDIR - Branch length generation


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

//=============================================================================
// branch length generation

float genBranch(float generate, float mean, float sdev)
{
    float blen = 0;
    const float maxlen = mean + sdev * 3.0;
    
    assert(mean != 0.0 || sdev != 0.0);

    while (blen <= 0 || blen > maxlen) {
        blen = normalvariate(mean, sdev);
	if (isnan(blen)) {
	    blen = mean;
	    break;
	}
    }
    
    return generate * blen;
}



// NOTE: currently generates zero branch length for freebranches...
void genSubtree(Tree *tree, Node *root,
                SpeciesTree *stree,
                int *ptree, int *pstree,
                int *recon, int *events, SpidirParams *params,
                float generate,
                ReconParams *reconparams)
{
    const int sroot = stree->root->name;

    
    if (events[root->name] != EVENT_DUP) {
        // single branch
                
        if (recon[root->name] != sroot) {
            // no midpoints
            reconBranch(root->name, tree, stree,
			ptree, pstree, recon, events, params,
                        reconparams);
            BranchParams bparam = getBranchParams(root->name, ptree, reconparams);            
            root->dist = genBranch(generate, bparam.alpha, bparam.beta);
        } else {
            // set branch length above sroot by exponential
            const float preratio = 0.05;
            root->dist = expovariate(1.0 / (generate * preratio));
        }
        
    } else {
        // multiple branches
                
        // set reconparams by traversing subtree
        ExtendArray<Node*> subnodes(0, tree->nnodes);
        getSubtree(root, events, &subnodes);
        
        ExtendArray<int> subnames(subnodes.size());
        for (int i=0; i<subnodes.size(); i++)
            subnames[i] = subnodes[i]->name;
        
        
        for (int i=0; i<subnodes.size(); i++) {
            reconBranch(subnodes[i]->name, tree, stree, ptree, pstree, 
                        recon, events, params, reconparams);
        }
        
        
        // propose a setting of midpoints
        reconparams->midpoints[root->name] = 1.0; // TODO: need to understand why this is here
        setRandomMidpoints(root->name, ptree, subnames, subnodes.size(),
                           recon, events, reconparams);
        
        // loop through all branches in subtree
        for (int j=0; j<subnodes.size(); j++) {
            Node *node = subnodes[j];

            if (recon[node->name] != sroot) {
                BranchParams bparam = getBranchParams(node->name, ptree, reconparams);
                node->dist = genBranch(generate, bparam.alpha, bparam.beta);
            } else {
                // set branch length above sroot by exponential
                const float preratio = 0.05;
		int i = 0;

		do {
		    node->dist = expovariate(1.0 / (generate * preratio));
		    i += 1;
		    if (i > 10) {
			node->dist = .001;
		    }
		} while (isnan(node->dist));
            }
        }
    }
}

// TODO: here is extra doc that should be worked in
// subnode == -1 && subnode == -1  --> generate whole tree
// subnode == -2 --> randomly pick node and child to resample

void generateBranchLengths(Tree *tree,
                           SpeciesTree *stree,
                           int *recon, int *events,
                           SpidirParams *params,
                           float generate, int subnode, int subchild)
{
    // generate a gene rate if it is requested (i.e. generate < 0)
    if (generate < 0.0)
        generate = gammavariate(params->gene_alpha, params->gene_beta);
    
    
    // determine reconciliation parameters
    ReconParams reconparams = ReconParams(tree->nnodes);
    determineFreeBranches(tree, stree, recon, events, generate,
                          &reconparams.unfold, 
                          &reconparams.unfolddist, 
                          reconparams.freebranches);
    
    
    // make array formats
    ExtendArray<int> ptree(tree->nnodes);
    ExtendArray<int> pstree(stree->nnodes);
    tree2ptree(tree, ptree);
    tree2ptree(stree, pstree);
    
    if (subnode != -1 && subchild != -1) {
        if (subnode == -2) {
            // pick random subnode
            do {
                subnode = irand(tree->nnodes);                
            } while (events[subnode] != EVENT_SPEC &&
		     subnode != tree->root->name);
            subchild = irand(2);
        }

        // regenerate only one subtree
        Node *node = (subnode == tree->root->name) ?
	    tree->root :
	    tree->nodes[subnode]->children[subchild];
        genSubtree(tree, node,
                   stree, ptree, pstree,
                   recon, events, params,
                   generate,
                   &reconparams);
    } else {    
        // loop through independent subtrees
        for (int i=0; i<tree->nnodes; i++) {
            if (events[i] == EVENT_SPEC || i == tree->root->name) {
                for (int j=0; j<2; j++) {
                    Node *node = tree->nodes[i]->children[j];
                    genSubtree(tree, node,
                               stree, ptree, pstree,
                               recon, events, params,
                               generate,
                               &reconparams);
                }
            }
        }
    }
}


void generateBranchLengths(int nnodes, int *ptree, 
                           int nsnodes, int *pstree,
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta,
                           float *dists)
{
    // create gene tree object
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);
    
    // create species tree object
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();  
    
    // build parameters object
    SpidirParams params(nsnodes, NULL, mu, sigma, alpha, beta);
    
    generateBranchLengths(&tree, &stree, recon, events, &params);
    
    // record distances into array
    tree.getDists(dists);
}
                         
} // namespace spidir
