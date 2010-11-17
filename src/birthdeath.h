#ifndef SPIDIR_BIRTHDEATH_H
#define SPIDIR_BIRTHDEATH_H

#include "Tree.h"
#include "spidir.h"

namespace spidir {

extern "C" {

int inumHistories(int ngenes);

double numHistories(int ngenes);

int inumTopologyHistories(Tree *tree);

double numTopologyHistories(Tree *tree);

void calcDoomTable(Tree *tree, float birth, float death, int maxdoom,
                   double *doomtable);

double birthDeathCount(int ngenes, float time, float birth, float death);
double birthDeathCounts(int start, int end, float time, 
                        float birth, float death);


double birthDeathTreePrior(Tree *tree, Tree *stree, int *recon, 
                          int *events, float birth, float death,
                          double *doomtable, int maxdoom);

double birthDeathTreePriorFull(Tree *tree, Tree *stree, int *recon, 
                              int *events, float birth, float death,
                              double *doomtable, int maxdoom);

double birthDeathTreeQuickPrior(Tree *tree, SpeciesTree *stree, int *recon, 
                                int *events, float birth, float death);


double sampleBirthWaitTime1(float T, float birth, float death);


}

} // namespace spidir

#endif // SPIDIR_BIRTHDEATH_H
