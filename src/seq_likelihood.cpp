/*=============================================================================

    SPIDIR    
    Maximum Likelihood Branch Length Estimation
    
    mldist.cpp
    started: Sun Jul  1 13:11:02 EDT 2007

=============================================================================*/

// c++ headers
#include <math.h>
#include <time.h>


// 3rd party
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multimin.h>

// spidir headers
#include "common.h"
#include "Matrix.h"
#include "seq_likelihood.h"
#include "parsimony.h"
#include "spidir.h"
#include "Tree.h"
#include "hky.h"



namespace spidir {

// prototype
template <class Model>
void calcLkTableRow(int seqlen, Model &model,
		    float *lktablea, float *lktableb, float *lktablec, 
		    float adist, float bdist);



//=============================================================================


template <class Model, class DModel>
class DistLikelihoodDeriv
{
public:
    DistLikelihoodDeriv(int seqlen, 
			Model *model, DModel *dmodel) :
        probs1(NULL),
        probs2(NULL),
        seqlen(seqlen),
        bgfreq(NULL),
        model(model),
	dmodel(dmodel)
    {
	probs3 = new float [seqlen*4];
	probs4 = new float [seqlen*4];
    }
    
    ~DistLikelihoodDeriv()
    {
	delete [] probs3;
	delete [] probs4;
    }

    void set_params(float *_probs1, float *_probs2, const float *_bgfreq)
    {
        probs1 = _probs1;
        probs2 = _probs2;
        bgfreq = _bgfreq;
    }

    double operator()(float t)
    {
        // trivial case
        if (t < 0)
            return INFINITY;
	

	// g(t, j)
	calcLkTableRow(seqlen, *model, probs1, probs2, probs3, t, 0);

	// g'(t, j)
	calcDerivLkTableRow(seqlen, *dmodel, probs1, probs2, probs4, t);


	// interate over sequence
	double dlogl = 0.0;
	for (int j=0; j<seqlen; j++) {
	    double sum1 = 0.0, sum2 = 0.0;
	    for (int k=0; k<4; k++) {
		sum1 += bgfreq[k] * probs3[matind(4,j,k)];
		sum2 += bgfreq[k] * probs4[matind(4,j,k)];
	    }
	    dlogl += sum2 / sum1;
	}
        
	return dlogl;
    }
    
    float *probs1;
    float *probs2;
    float *probs3;
    float *probs4;
    int seqlen;
    const float *bgfreq;
    Model *model;
    DModel *dmodel;
};



template <class Model, class DModel, class D2Model>
class DistLikelihoodDeriv2
{
public:
    DistLikelihoodDeriv2(int seqlen, 			
                         Model *model, DModel *dmodel, D2Model *d2model) :
        seqlen(seqlen),
        model(model),
	dmodel(dmodel),
        d2model(d2model)
    {
	probs3 = new float [seqlen*4];
	probs4 = new float [seqlen*4];
        probs5 = new float [seqlen*4];
    }
    
    ~DistLikelihoodDeriv2()
    {
	delete [] probs3;
	delete [] probs4;
        delete [] probs5;
    }

    void set_params(float *_probs1, float *_probs2, const float *_bgfreq)
    {
        probs1 = _probs1;
        probs2 = _probs2;
        bgfreq = _bgfreq;
    }


    double operator()(float t)
    {
        // trivial case
        if (t < 0)
            return -INFINITY;
	

	// g(t, j)
	calcLkTableRow(seqlen, *model, probs1, probs2, probs3, t, 0);

	// g'(t, j)
	calcDerivLkTableRow(seqlen, *dmodel, probs1, probs2, probs4, t);

	// g''(t, j)
	calcDerivLkTableRow(seqlen, *d2model, probs1, probs2, probs5, t);


	// interate over sequence
	double d2logl = 0.0;
	for (int j=0; j<seqlen; j++) {
	    double g = 0.0, dg = 0.0, d2g = 0.0;
	    for (int k=0; k<4; k++) {
		g += bgfreq[k] * probs3[matind(4,j,k)];
		dg += bgfreq[k] * probs4[matind(4,j,k)];
                d2g += bgfreq[k] * probs5[matind(4,j,k)];
	    }
	    d2logl += - dg*dg/(g*g) + d2g/g;
	}
        
	return d2logl;
    }
    
    float *probs1;
    float *probs2;
    float *probs3;
    float *probs4;
    float *probs5;
    int seqlen;
    const float *bgfreq;
    Model *model;
    DModel *dmodel;
    D2Model *d2model;
};


// Find the maximum likelihood estimate (MLE) of the distance (time) between 
// two sequences (represented probabilistically as probs1 and probs2)
//
//  bgfreq = background frequency of bases
//  model = sequence evolution model
//
template <class Model>
float mleDistance(float *probs1, float *probs2, int seqlen, 
                  const float *bgfreq, Model &model, 
                  float t0=.001, float t1=1, float step=.0001)
{
    typename Model::Deriv *dmodel = model.deriv();

    DistLikelihoodDeriv<Model, typename Model::Deriv> df(
        probs1, probs2, seqlen, bgfreq, &model, dmodel);
    float mle = bisectRoot(df, t0, t1, .0001);

    delete dmodel;
    return mle;
}

extern "C" {

/*

  lk = prod_j sum_k bgfreq[k] * lktable[root][j,k] 
     = prod_j sum_k bgfreq[k] * (sum_x P(x|k, t_a) lktable[a][j,x]) *
                                (sum_y P(y|k, t_b) lktable[b][j,y])

 */
float branchLikelihoodHky(float *probs1, float *probs2, int seqlen, 
			  const float *bgfreq, float kappa, float t)
{

    HkyModel hky(bgfreq, kappa);

    // allocate precomputed probability terms
    float *probs3 = new float [seqlen*4];
    
    calcLkTableRow(seqlen, hky, probs1, probs2, probs3, 0, t);

    float logl = 0.0;

    // interate over sequence
    for (int j=0; j<seqlen; j++) {
	float sum = 0.0;
	for (int k=0; k<4; k++)
	    sum += bgfreq[k] * probs3[matind(4,j,k)];
	logl += logf(sum);
    }

    // free probability table
    delete [] probs3;

    return logl;
}


float branchLikelihoodHkyDeriv(float *probs1, float *probs2, int seqlen, 
                               const float *bgfreq, float kappa, float t)
{
    HkyModel hky(bgfreq, kappa);
    HkyModelDeriv dhky(bgfreq, kappa);
    DistLikelihoodDeriv<HkyModel, HkyModelDeriv> df(seqlen, &hky, &dhky);
    df.set_params(probs1, probs2, bgfreq);
    return df(t);
}

float branchLikelihoodHkyDeriv2(float *probs1, float *probs2, int seqlen, 
                               const float *bgfreq, float kappa, float t)
{
    HkyModel hky(bgfreq, kappa);
    HkyModelDeriv dhky(bgfreq, kappa);
    HkyModelDeriv2 d2hky(bgfreq, kappa);
    DistLikelihoodDeriv2<HkyModel, HkyModelDeriv, HkyModelDeriv2> d2f
	(seqlen, &hky, &dhky, &d2hky);
    d2f.set_params(probs1, probs2, bgfreq);
    return d2f(t);
}

float mleDistanceHky(float *probs1, float *probs2, int seqlen, 
		     const float *bgfreq, float kappa,
		     float t0, float t1)
{
    HkyModel hky(bgfreq, kappa);
    HkyModelDeriv dhky(bgfreq, kappa);
    DistLikelihoodDeriv<HkyModel, HkyModelDeriv> df(seqlen, &hky, &dhky);
    df.set_params(probs1, probs2, bgfreq);
    
    return bisectRoot(df, t0, t1, .0001);
}

} // extern "C"


//=============================================================================
// Dynamic programming of conditional likelihood


class LikelihoodTable 
{
public:

    LikelihoodTable(int nnodes, int seqlen) :
        nnodes(nnodes),
        seqlen(seqlen)
    {
        // allocate conditional likelihood dynamic programming table
        lktable = new float* [nnodes];
        for (int i=0; i<nnodes; i++)
            lktable[i] = new float [4 * seqlen];
    }

    ~LikelihoodTable()
    {
        // cleanup
        for (int i=0; i<nnodes; i++)
            delete [] lktable[i];
    }

    float **lktable;

    int nnodes;
    int seqlen;
    
};



/*

    From: Inferring Phylogenies. Felsenstein. p 254.
    
    Baisc recursion of conditional likelihoods
        lktable[c][j,k] = (sum_x P(x|k, t_a) lktable[a][j,x]) *
                          (sum_y P(y|k, t_b) lktable[b][j,y])
    
    where c is a node, with child nodes a and b
          t_a, t_b are the branch lengths of branches a and b
          j indexes sites
          k,x,y are DNA bases (e.g. A=1, C=2, G=3, T=4)

*/


// conditional likelihood recurrence
template <class Model>
void calcLkTableRow(int seqlen, Model &model,
		    float *lktablea, float *lktableb, float *lktablec, 
		    float adist, float bdist)
{
    float atransmat[16];
    float btransmat[16];
    
    // build transition matrices
    model.getMatrix(adist, atransmat);
    model.getMatrix(bdist, btransmat);
    
    // iterate over sites
    for (int j=0; j<seqlen; j++) {
        const float *terma = &lktablea[matind(4, j, 0)];
        const float *termb = &lktableb[matind(4, j, 0)];
        
        for (int k=0; k<4; k++) {
            const float *aptr = &atransmat[4*k];
            const float *bptr = &btransmat[4*k];
            
            // sum_x P(x|k, t_a) lktable[a][j,x]
            const float prob1 = aptr[0] * terma[0] +
                                aptr[1] * terma[1] +
                                aptr[2] * terma[2] +
                                aptr[3] * terma[3];

            // sum_y P(y|k, t_b) lktable[b][j,y]
            const float prob2 = bptr[0] * termb[0] +
                                bptr[1] * termb[1] +
                                bptr[2] * termb[2] +
                                bptr[3] * termb[3];
            

            lktablec[matind(4, j, k)] = prob1 * prob2;
        }
    }
}


// conditional likelihood recurrence
template <class DModel>
void calcDerivLkTableRow(int seqlen, DModel &dmodel,
			 float *lktablea, float *lktableb, float *lktablec, 
			 float adist)
{
    float btransmat[16];
    
    // build transition matrix
    dmodel.getMatrix(adist, btransmat);
    
    // iterate over sites
    for (int j=0; j<seqlen; j++) {
        const float *terma = &lktablea[matind(4, j, 0)];
        const float *termb = &lktableb[matind(4, j, 0)];
        
        for (int k=0; k<4; k++) {
            //float *aptr = &atransmat[4*k];
            const float *bptr = &btransmat[4*k];
            
            // sum_y P(y|k, t_b) lktable[b][j,y]
            const float prob2 = bptr[0] * termb[0] +
                                bptr[1] * termb[1] +
                                bptr[2] * termb[2] +
                                bptr[3] * termb[3];
            
            lktablec[matind(4, j, k)] = terma[k] * prob2;
        }
    }
}


// initialize the condition likelihood table
template <class Model>
void calcLkTable(float** lktable, Tree *tree, 
                 int nseqs, int seqlen, char **seqs, Model &model)
{
    // recursively calculate cond. lk. of internal nodes
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);
    
    for (int l=0; l<nodes.size(); l++) {
        Node *node = nodes[l];
        int i = node->name;
        
        if (node->isLeaf()) {
            // initialize leaves from sequence
        
            // iterate over sites
            for (int j=0; j<seqlen; j++) {
                int base = dna2int[int(seqs[i][j])];

                if (base == -1) {
                    // handle gaps
                    lktable[i][matind(4, j, 0)] = 1.0;
                    lktable[i][matind(4, j, 1)] = 1.0;
                    lktable[i][matind(4, j, 2)] = 1.0;
                    lktable[i][matind(4, j, 3)] = 1.0;
                } else {
                    lktable[i][matind(4, j, 0)] = 0.0;
                    lktable[i][matind(4, j, 1)] = 0.0;
                    lktable[i][matind(4, j, 2)] = 0.0;
                    lktable[i][matind(4, j, 3)] = 0.0;

                    lktable[i][matind(4, j, base)] = 1.0;
                }
            }
        } else {
            // compute internal nodes from children
            Node *node1 = node->children[0];
            Node *node2 = node->children[1];
            
            calcLkTableRow(seqlen, model, 
                           lktable[node1->name], 
			   lktable[node2->name], 
			   lktable[node->name],
                           node1->dist, node2->dist);
        }
    }
}


// calculate log(P(D | T, B))
template <class Model>
float getTotalLikelihood(float** lktable, Tree *tree, 
                         int seqlen, Model &model, const float *bgfreq)
{
    // integrate over the background base frequency
    const float *rootseq = lktable[tree->root->name];
    float lk = 0.0;
    for (int k=0; k<seqlen; k++) {
        float prob = 0.0;
        for (int x=0; x<4; x++)
            prob += bgfreq[x] * rootseq[matind(4, k, x)];
        lk += logf(prob);
    }

    // return log likelihood
    return lk;
}



template <class Model>
float calcSeqProb(Tree *tree, int nseqs, char **seqs, 
                  const float *bgfreq, Model &model)
{
    int seqlen = strlen(seqs[0]);
    
    LikelihoodTable table(tree->nnodes, seqlen);
    calcLkTable(table.lktable, tree, nseqs, seqlen, seqs, model);
    float logl = getTotalLikelihood(table.lktable, tree, seqlen, model, bgfreq);
    
    return logl;
}

extern "C" {

float calcSeqProbHky(Tree *tree, int nseqs, char **seqs, 
                     const float *bgfreq, float ratio)
{
    HkyModel hky(bgfreq, ratio);
    return calcSeqProb(tree, nseqs, seqs, bgfreq, hky);
}

} // extern "C"


//=============================================================================
// find MLE branch lengths


// TODO: need to think about more carefully
void getRootOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL)
{
    if (!node) {
        // start at tree root (but don't include root)
        getRootOrder(tree, nodes, tree->root->children[0]);
        getRootOrder(tree, nodes, tree->root->children[1]);
    } else {
        // record pre-process
        if (node->parent == tree->root) {
            nodes->append(tree->root->children[0]);
            nodes->append(tree->root->children[1]);
        } else {
            nodes->append(node);
            nodes->append(node->parent);
        }

        // recurse
        for (int i=0; i<node->nchildren; i++)
            getRootOrder(tree, nodes, node->children[i]);
        
        /*
        // record post-process
        if (!node->isLeaf()) {
            if (node->parent == tree->root) {
                nodes->append(tree->root->children[0]);
                nodes->append(tree->root->children[1]);
            } else {
                nodes->append(node);
                nodes->append(node->parent);
            }
        }*/
    }
}

template <class Model>
class MLBranchAlgorithm
{
public:

    MLBranchAlgorithm(Tree *tree, int seqlen, Model *model) :
	table(tree->nnodes, seqlen),
        model(model),
        dmodel(model->deriv()),
        d2model(dmodel->deriv()),
        lk_deriv(seqlen, model, dmodel),
	lk_deriv2(seqlen, model, dmodel, d2model)
    {

        // setup branch length optimizer
        opt = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);

        opt_func.f = &branch_f;
        opt_func.df = &branch_df;
        opt_func.fdf = &branch_fdf;
        opt_func.params = this;
    }

    ~MLBranchAlgorithm()
    {
        gsl_root_fdfsolver_free(opt);
        delete dmodel;
    }


    static double branch_f(double x, void *params)
    {       
        double y;
        if (x < minx || x > maxx)
            y = ((MLBranchAlgorithm*) params)->lk_deriv(minx);
        else
            y = ((MLBranchAlgorithm*) params)->lk_deriv(x);
        //printf("x %f y %f\n", x, y);
        return y;
    }

    static double branch_df(double x, void *params)
    {
        double dy;
        if (x < minx || x > maxx)
            dy = ((MLBranchAlgorithm*) params)->lk_deriv(minx) / (x - minx);
        else
            dy = ((MLBranchAlgorithm*) params)->lk_deriv2(x);
        
        //printf("x %f dy %f\n", x, dy);
        return dy;
    }

    static void branch_fdf(double x, void *params, 
                           double *f, double *df)
    {
        *f = branch_f(x, params); //((MLBranchAlgorithm*) params)->lk_deriv(x);
        *df = branch_df(x, params); //((MLBranchAlgorithm*) params)->lk_deriv2(x);
        //printf("x %f y %f dy %f\n", x, *f, *df);
    }


    // old branch fitting
    float fitBranch2(Tree *tree, const float *bgfreq, float initdist)
    {
        Node *node1 = tree->root->children[0];
        Node *node2 = tree->root->children[1];

        lk_deriv.set_params(table.lktable[node1->name], 
                            table.lktable[node2->name], 
                            bgfreq);

        return bisectRoot(lk_deriv, max(initdist*0.0, 0.0), 
                          max(initdist*10.0, 0.00001), .0001);

    }

    float fitBranch(Tree *tree, const float *bgfreq, float initdist)
    {
        double r = initdist, r0 = initdist;
        const double esp = 1e-3;

        if (initdist < 0)
            initdist = .0001;

        //printf("initdist %f\n", initdist);

        Node *node1 = tree->root->children[0];
        Node *node2 = tree->root->children[1];
        
        lk_deriv.set_params(table.lktable[node1->name], 
                            table.lktable[node2->name], 
                            bgfreq);
        
        lk_deriv2.set_params(table.lktable[node1->name], 
                             table.lktable[node2->name], 
                             bgfreq);

        gsl_root_fdfsolver_set(opt, &opt_func, initdist);

        int status = GSL_CONTINUE;
        const int maxiter = 10;
        int iter;
        for (iter=0; iter<maxiter && status==GSL_CONTINUE; iter++) {
            //printf("root %f\n", r);

            // do one iteration
            status = gsl_root_fdfsolver_iterate(opt);
            r0 = r;
            r = gsl_root_fdfsolver_root(opt);
            if (status)
                break;
            status = gsl_root_test_delta(r, r0, 0, esp);
        }
        //printf("root done iters=%d r=%f\n", iter, r);

        if (iter == maxiter || status != GSL_SUCCESS) {
            r = fitBranch2(tree, bgfreq, initdist);
        }

        //printf("bisect: %f\n\n", r);

        // return final branch length
        return r;
    }
        



    float fitBranches(Tree *tree, int nseqs, int seqlen, char **seqs, 
		      const float *bgfreq, 
		      ExtendArray<Node*> &rootingOrder)
    {
	float logl = -INFINITY;
	float **lktable = table.lktable;


	for (int i=0; i<rootingOrder.size(); i+=2) {
	    // remembering old children of root
	    Node *oldnode1 = tree->root->children[0];
	    Node *oldnode2 = tree->root->children[1];

	    // choose new root
	    tree->reroot(rootingOrder[i], rootingOrder[i+1]);
	    if (tree->root->children[0] == oldnode1 &&
		tree->root->children[1] == oldnode2)
		continue;	    

	    // rebuild invaild likelihood values
	    // determine starting node to rebuild
	    Node *ptr;
	    if (oldnode1->parent == oldnode2)
		ptr = oldnode1;
	    else if (oldnode2->parent == oldnode1)
		ptr = oldnode2;
	    else
		assert(0);

	    // walk up to root of tree, rebuilding conditional likelihoods
	    for (; ptr; ptr = ptr->parent) {
		if (!ptr->isLeaf())
		    calcLkTableRow(seqlen, *model, 
				   lktable[ptr->children[0]->name], 
				   lktable[ptr->children[1]->name], 
				   lktable[ptr->name],
				   ptr->children[0]->dist, 
				   ptr->children[1]->dist);
	    }

	    // get total probability before branch length change
	    float loglBefore = getTotalLikelihood(lktable, tree, 
                                                  seqlen, *model, bgfreq);

            Node *node1 = tree->root->children[0];
            Node *node2 = tree->root->children[1];
            float initdist = node1->dist + node2->dist;

            // find new MLE branch length for root branch
            float mle = fitBranch(tree, bgfreq, initdist);
	    node1->dist = mle / 2.0;
	    node2->dist = mle / 2.0;
	

	    // recompute the root node row in lktable
	    calcLkTableRow(seqlen, *model, 
			   lktable[node1->name], 
			   lktable[node2->name], 
			   lktable[tree->root->name],
			   node1->dist, 
			   node2->dist);
	
	    // get total probability after branch change    
	    logl = getTotalLikelihood(lktable, tree, 
				      seqlen, *model, bgfreq);
	
	    // don't accept a new branch length if it lowers total likelihood
	    if (logl < loglBefore) {
	         // revert
	    	 node1->dist = initdist / 2.0;
	    	 node2->dist = initdist / 2.0;
	    	 logl = loglBefore;
	    	 //assert(0);
	    }
        
	    printLog(LOG_HIGH, "hky: lk=%f\n", logl);        
	}

	return logl;
    }

    const static double minx = 0.0000001;
    const static double maxx = 10.0;

    gsl_root_fdfsolver *opt;
    gsl_function_fdf opt_func;

    LikelihoodTable table;
    Model *model;
    typename Model::Deriv *dmodel;
    typename Model::Deriv::Deriv *d2model;
    
    DistLikelihoodDeriv<Model, typename Model::Deriv> lk_deriv;
    DistLikelihoodDeriv2<Model, typename Model::Deriv, 
                         typename Model::Deriv::Deriv> lk_deriv2;
};



// NOTE: assumes binary Tree
template <class Model>
float findMLBranchLengths(Tree *tree, int nseqs, char **seqs, 
                          const float *bgfreq, Model &model,
                          int maxiter=10)
{
    // timing
    Timer timer;
    

    int seqlen = strlen(seqs[0]);
    MLBranchAlgorithm<Model> mlalg(tree, seqlen, &model);

    float lastLogl = -INFINITY, logl = -INFINITY;    
    const float converge = logf(1.02);
    
    // initialize the condition likelihood table
    calcLkTable(mlalg.table.lktable, tree, nseqs, seqlen, seqs, model);
    
    
    // determine rooting order
    ExtendArray<Node*> rootingOrder(0, 2*tree->nnodes);
    getRootOrder(tree, &rootingOrder);
    
    // remember original rooting for restoring later
    Node *origroot1 = tree->root->children[0];
    Node *origroot2 = tree->root->children[1];
    
    // iterate over branches improving each likelihood
    for (int j=0; j<maxiter; j++) {
        printLog(LOG_HIGH, "hky: iter %d\n", j);  

	logl = mlalg.fitBranches(tree, nseqs, seqlen, seqs, bgfreq, 
				 rootingOrder);
        
        // determine whether logl has converged
        float diff = fabs(logl - lastLogl);
        if (diff < converge) {
            printLog(LOG_HIGH, "hky: diff = %f < %f\n", diff, converge);
            break;
        } else {
            printLog(LOG_HIGH, "hky: diff = %f > %f\n", diff, converge);
        }
        lastLogl = logl;
    }
    
    // restore original rooting
    tree->reroot(origroot1, origroot2);

    
    printLog(LOG_MEDIUM, "mldist time: %f\n",  timer.time());
    
    return logl;
}


float findMLBranchLengthsHky(Tree *tree, int nseqs, char **seqs, 
                             const float *bgfreq, float ratio, int maxiter)
{
    HkyModel hky(bgfreq, ratio);
    return findMLBranchLengths(tree, nseqs, seqs, bgfreq, hky, maxiter);
}



extern "C" {

float findMLBranchLengthsHky(int nnodes, int *ptree, int nseqs, char **seqs, 
                             float *dists, const float *bgfreq, float ratio, 
                             int maxiter, bool parsinit)
{
    //int seqlen = strlen(seqs[0]);
        
    // create tree objects
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);

    if (parsinit)
        parsimony(&tree, nseqs, seqs);
    
    float logl = findMLBranchLengthsHky(&tree, nseqs, seqs, bgfreq, 
                                        ratio, maxiter);
    tree.getDists(dists);
    
    return logl;
}


} // extern C


} // namespace spidir
