/* This file has code for solving for stochastic user equilibrium using
 * the method of successive averages with fixed step size.
 */

#include "convexcombination.h"

#define MAX_TIME 3600 /* Maximum run time, in seconds */
#define MAX_ITERATIONS 100 /* Maximum # of iterations */
#define LINK_FLOW_TOLERANCE 1e-3 /* Stop if average link flow is this close
                                    to target */

/* Main function for method of successive averages with fixed
 * step size lambda
 */
void SUE_MSA(network_type *network, double theta, double lambda) {
    bool converged = FALSE;
    int iteration = 0;
    double elapsedTime = 0, diff = INFINITY;
    long numBushLinks, numPaths;
    bushes_type *bushes = NULL;
    declareVector(double, target, network->numArcs);
    clock_t stopTime = clock(); /* used for timing */

    initializeSolution(network, &bushes, theta, &numBushLinks, &numPaths);
    elapsedTime += ((double)(clock() - stopTime)) / CLOCKS_PER_SEC;
    displayMessage(MEDIUM_NOTIFICATIONS, "%ld bush links, %ld paths\n",
                   numBushLinks, numPaths);
    displayMessage(LOW_NOTIFICATIONS, "Initialization done in %.3f s.\n",
                   elapsedTime);
    stopTime = clock();
    while (converged == FALSE) {
        updateLinkCosts(network);
        calculateTarget(network, bushes, target, theta);
        diff = avgFlowDiff(network, target);
        elapsedTime += ((double)(clock() - stopTime)) / CLOCKS_PER_SEC;
        displayMessage(LOW_NOTIFICATIONS, "Iteration %d: "
                                          "flow diff %.3f, "
                                          "time %.3f\n",
                                          iteration,
                                          diff,
                                          elapsedTime);
        stopTime = clock();
        if (elapsedTime > MAX_TIME) converged = TRUE;
        if (iteration >= MAX_ITERATIONS) converged = TRUE;
        if (diff < LINK_FLOW_TOLERANCE) converged = TRUE;
        if (converged == TRUE) break;

        shiftFlows(network, target, lambda);
        iteration++;
    }
    deleteVector(target);
    deleteBushes(bushes);
}

/* 
 * Adjust link flows by taking a step of a given size in the direction 
 * of a given target.
 */
void shiftFlows(network_type *network, double *target, double stepSize) {
    int ij;
    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].flow += stepSize * (target[ij] 
                                              - network->arcs[ij].flow);
    }
}

/*
 * Compute target link flows by using Dial's method for each origin,
 * then summing the flows into a single array.
 */
void calculateTarget(network_type *network, bushes_type *bushes,
                     double *target, double theta) {
    int r, ij;
    for (ij = 0; ij < network->numArcs; ij++) {
        target[ij] = 0;
    }
    for (r = 0; r < network->numZones; r++) {
        dialFlows(network, bushes, r, theta);
        for (ij = 0; ij < network->numArcs; ij++) {
            target[ij] += bushes->flow[ij];
        }
    }
}

/*
 * Compute the average absolute difference between the current flow
 * vector and the target vector (to be used in a convergence check).
 */
double avgFlowDiff(network_type *network, double *target) {
    int ij;
    double total = 0;
    for (ij = 0; ij < network->numArcs; ij++) {
        total += fabs(network->arcs[ij].flow - target[ij]);
    }
    return total / network->numArcs;
}

/*
 * Generate an initial feasible solution and set up the bush data structures
 * for Dial's method in subsequent iterations.  Also compute the number
 * of bush links and number of bush paths to see how much space/calculation
 * is saved by using Dial's method rather than directly using the logit
 * formula.
 */
void initializeSolution(network_type *network, bushes_type **bushes,
                          double theta, long *numBushLinks, long *numPaths) {

    int r, ij;
    declareVector(double, target, network->numArcs);
    *bushes = initializeBushes(network);
    *numBushLinks = 0;
    *numPaths = 0;
    for (r = 0; r < network->numZones; r++) {
        *numBushLinks += (*bushes)->numBushLinks[r];
        *numPaths += (*bushes)->numBushPaths[r];
    }

    /* Compute initial solution */
    calculateTarget(network, *bushes, target, theta);
    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].flow = target[ij];
        //printf("%d %d has %f\n", network->arcs[ij].tail+1, network->arcs[ij].head+1, network->arcs[ij].flow);
    }
    deleteVector(target);
}
