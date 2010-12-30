#ifndef SPIDIR_TOP_PRIOR_H
#define SPIDIR_TOP_PRIOR_H

#include "Tree.h"
#include "spidir.h"

namespace spidir {

extern "C" {

void calcDoomTable(Tree *tree, float birth, float death, int maxdoom,
                   double *doomtable);

void getSpecSubtree(Node *node, Node *snode, int *recon, int *events,
                    ExtendArray<Node*> &nodes);

double birthDeathTreePrior(Tree *tree, Tree *stree, int *recon, 
                          int *events, float birth, float death,
                          double *doomtable, int maxdoom);

double birthDeathTreePriorFull(Tree *tree, Tree *stree, int *recon, 
                              int *events, float birth, float death,
                              double *doomtable, int maxdoom);

double birthDeathTreeQuickPrior(Tree *tree, SpeciesTree *stree, int *recon, 
                                int *events, float birth, float death);




}

} // namespace spidir

#endif // SPIDIR_TOP_PRIOR_H
