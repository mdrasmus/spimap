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
                   float *doomtable);

float birthDeathCount(int ngenes, float time, float birth, float death);
float birthDeathCounts(int start, int end, float time, 
                       float birth, float death);

double birthDeathTreeCounts(Tree *tree, int nspecies, int *counts, 
                            float birth, float death, int maxgene,
                            int rootgene, double **tab=NULL);

double birthDeathForestCounts(Tree *tree, int nspecies, int nfams,
                              int **counts, int *mult,
                              float birth, float death, int maxgene,
                              int rootgene, double **tab=NULL);

float birthDeathDensity(float *times, int ntimes, float maxtime, 
                        float birth, float death);

float birthDeathTreePrior(Tree *tree, Tree *stree, int *recon, 
                          int *events, float birth, float death,
                          float *doomtable, int maxdoom);

float birthDeathTreePriorFull(Tree *tree, Tree *stree, int *recon, 
                              int *events, float birth, float death,
                              float *doomtable, int maxdoom);

float birthDeathTreeQuickPrior(Tree *tree, SpeciesTree *stree, int *recon, 
                               int *events, float birth, float death);

void sampleDupTimes(Tree *tree, Tree *stree, int *recon, int *events,
                    float birth, float death);


float sampleBirthWaitTime1(float T, float birth, float death);


}

} // namespace spidir

#endif // SPIDIR_BIRTHDEATH_H
