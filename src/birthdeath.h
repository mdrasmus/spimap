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

void calcDoomTable(Tree *tree, float birthRate, float deathRate, int maxdoom,
                   float *doomtable);

float birthDeathCount(int ngenes, float time, float birthRate, float deathRate);

float birthDeathDensity(float *times, int ntimes, float maxtime, 
                        float birthRate, float deathRate);

float birthDeathTreePrior(Tree *tree, Tree *stree, int *recon, 
                          int *events, float birthRate, float deathRate,
                          float *doomtable, int maxdoom);

float birthDeathTreePriorFull(Tree *tree, Tree *stree, int *recon, 
                              int *events, float birthRate, float deathRate,
                              float *doomtable, int maxdoom);

float birthDeathTreeQuickPrior(Tree *tree, SpeciesTree *stree, int *recon, 
                               int *events, float birthRate, float deathRate);

void sampleDupTimes(Tree *tree, SpeciesTree *stree, int *recon, int *events,
                    float birthRate, float deathRate);


}

} // namespace spidir

#endif // SPIDIR_BIRTHDEATH_H
