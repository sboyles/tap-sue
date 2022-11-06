#ifndef CONVEXCOMBINATIONS_H
#define CONVEXCOMBINATIONS_H

#include <limits.h>
#include <stdio.h>
#include <time.h>
#include "fileio.h"
#include "bush.h"
#include "networks.h"
#include "utils.h"

void SUE_MSA(network_type *network, double theta, double lambda);
void shiftFlows(network_type *network, double *target, double stepSize);
void calculateTarget(network_type *network, bushes_type *bushes,
                       double *target, double theta);
double avgFlowDiff(network_type *network, double *target);
void initializeSolution(network_type *network, bushes_type *bushes,
                          double theta, long *numBushLinks, long *numPaths);
#endif
