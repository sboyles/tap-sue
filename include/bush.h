/*
 * bush.h -- This is the header file for the bush data structure.
 * These structures are the key for handling Dial's method efficiently.
 *
 * Each origin is associated with a "bush" of "reasonable links", indicated
 * by origin-specific forward and reverse star lists, and an origin-specific
 * topological order.
 *
 * This file also includes functions for implementing Dial's method.
 */

#ifndef BUSH_H
#define BUSH_H

#include <math.h>
#include "networks.h"
#include "datastructures.h"
#include "utils.h"

#define MIN_LINK_COST 1e-6 /* Ensure links have strictly positive cost
                              for finding initial bushes */

/*
 * bushes_type: Stores the data for ALL bushes associated with the network
 *
 * The following members are *shared* across all bushes, and are associated
 * with whatever bush is currently being operated on.  To save memory, this
 * information is overwritten when we move to another bush.
 *  SPcost -- array of shortest path costs, indexed by node ID
 *  flow -- array of bush flows, indexed by link ID.
 *  nodeFlow -- array of total flow through each node in the bush, indexed by
 *              node ID.
 *  weight -- array of bush link weights, indexed by link ID.
 *  nodeWeight -- array of total weight at each node, indexed by node ID.
 *  likelihood -- array of link likelihoods, indexed by link ID.
 *
 * The following members are stored *separately* for each bush, and contain
 * information about bushes which is persistent even when other bushes are
 * being operated on.
 *  bushOrder -- stores the inverse topological order for each bush, that is,
 *               bushOrder[origin][index] gives the link ID of the node which
 *               is in the index-th position topologically, for bush 'origin'.
 *  bushForwardStar -- an arcList storing the reasonable links leaving a
 *                     particular node, indexed by origin and node.
 *  bushReverseStar -- like bushForwardStar, but for entering reasonable links
 *
 * The following are statistics about the bushes themselves:
 *  numBushLinks -- the number of reasonable links for a given origin
 *  numBushPaths -- the number of reasonable paths for a given origin
 */

typedef struct bushes_type {
    double *SPcost; /* [node] */
    double *flow; /* [link] */
    double *nodeFlow; /* [node] */
    double *weight; /* [link] */
    double *nodeWeight; /* [node] */
    double *likelihood; /* [link] */
    int **bushOrder; /* [origin][node] */
    arcList **bushForwardStar; /* [origin][node] */
    arcList **bushReverseStar; /* [origin][node] */
    network_type *network; /* Points back to the corresponding network */
    long *numBushLinks; /* [origin] */
    unsigned long long int *numBushPaths; /* [origin] */
} bushes_type;

bushes_type *initializeBushes(network_type *network);
void deleteBushes(bushes_type *bushes);

void bushTopologicalOrder(int origin, network_type *network,
                          bushes_type *bushes);
void bushShortestPath(network_type *network, bushes_type *bushes, int origin);
void dialFlows(network_type *network, bushes_type *bushes, int origin,
               double theta);
#endif
