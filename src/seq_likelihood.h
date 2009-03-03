#ifndef SPIDIR_SEQ_LIKELIHOOD_H
#define SPIDIR_SEQ_LIKELIHOOD_H


#include "Tree.h"

namespace spidir {


float findMLBranchLengthsHky(Tree *tree, int nseqs, char **seqs, 
                            const float *bgfreq, float ratio, 
                            int maxiter=100, int samples=0);

template <class Model>
float getTotalLikelihood(ExtendArray<float*> &lktable, Tree *tree, 
                         int nseqs, int seqlen, char **seqs, Model &model,
                         const float *bgfreq);

extern "C" {

void makeHkyMatrix(const float *bgfreq, float ratio, float t, float *matrix);

float branchLikelihoodHky(float *probs1, float *probs2, int seqlen, 
			  const float *bgfreq, float kappa, float t);

float deriveBranchLikelihoodHky(float *probs1, float *probs2, int seqlen, 
				const float *bgfreq, float kappa, float t);

float mleDistanceHky(float *probs1, float *probs2, int seqlen, 
		     const float *bgfreq, float ratio,
		     float t0, float t1, float step);

float calcSeqProbHky(Tree *tree, int nseqs, char **seqs, 
                     const float *bgfreq, float ratio);

float findMLBranchLengthsHky(int nnodes, int *ptree, int nseqs, char **seqs, 
                             float *dists, const float *bgfreq, float ratio, 
                             int maxiter, bool parsinit=false);

} // extern "C"

} // namespace spidir

#endif // SPIDIR_SEQ_LIKELIHOOD_H
