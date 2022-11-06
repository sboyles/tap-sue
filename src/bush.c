/*
 * bush.c -- Contains all of the components of Algorithm B, including bush
 * creation, updating, and flow shifting.
 *
 * See comments on bush.h for descriptions of the data structures used for
 * bushes.
 */
#include "bush.h"

/* Initialize bushes based on free-flow travel times. */
bushes_type *initializeBushes(network_type *network) {
    int r, curnode, i, j, ij;
    arcListElt *curArc;
    bushes_type *bushes = newScalar(bushes_type);

    bushes->SPcost = newVector(network->numNodes, double);
    bushes->flow = newVector(network->numArcs, double);
    bushes->nodeFlow = newVector(network->numNodes, double);
    bushes->weight = newVector(network->numArcs, double);
    bushes->nodeWeight = newVector(network->numNodes, double);
    bushes->likelihood = newVector(network->numArcs, double);
    
    bushes->bushForwardStar
            = newMatrix(network->numZones, network->numNodes, arcList);
    bushes->bushReverseStar
            = newMatrix(network->numZones, network->numNodes, arcList);

    bushes->numBushLinks = newVector(network->numZones, long);
    bushes->numBushPaths = newVector(network->numZones, long);

    declareVector(long, pathCount, network->numNodes);

    /* Identify reasonable links based on distance from origin using
     * free-flow costs.  Ensure these are strictly positive to prevent
     * issues with zero-cost links. */
    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].cost = max(MIN_LINK_COST,
                                     network->arcs[ij].freeFlowTime
                                         + network->arcs[ij].fixedCost);
    }
    for (r = 0; r < network->numZones; r++) {
        bushes->numBushLinks[r] = 0;
        for (i = 0; i < network->numNodes; i++) {
            initializeArcList(&(bushes->bushForwardStar[r][i]));
            initializeArcList(&(bushes->bushReverseStar[r][i]));
        }

        /* Identify reasonable links */
        shortestPath(r, bushes->SPcost, network);
        for (ij = 0; ij < network->numArcs; ij++) {
            i = network->arcs[ij].tail;
            j = network->arcs[ij].head;
            if (bushes->SPcost[i] < bushes->SPcost[j]) {
                bushes->numBushLinks[r]++;
                insertArcList(&(bushes->bushForwardStar[r][i]),
                              &(network->arcs[ij]),
                              bushes->bushForwardStar[r][i].tail);
                insertArcList(&(bushes->bushReverseStar[r][j]),
                              &(network->arcs[ij]),
                              bushes->bushForwardStar[r][j].tail);
            }
        }
        bushTopologicalOrder(r, network, bushes);

        /* Compute number of paths in the bush */
        bushes->numBushPaths[r] = 0;
        pathCount[r] = 1;
        for (curnode = 1; curnode < network->numNodes; curnode++) {
            j = bushes->bushOrder[r][curnode];
            pathCount[j] = 0;
            for (curArc = bushes->bushReverseStar[r][i].head;
                 curArc != NULL;
                 curArc = curArc->next)
            {
                i = curArc->arc->tail;
                pathCount[j] += pathCount[i];
            }
            if (j < network->numZones) {
                bushes->numBushPaths[r] += pathCount[j];
            }
        }

    }

    bushes->network = network;
    deleteVector(pathCount);
    return bushes;
}

/* Free memory associated with bush set. */
void deleteBushes(bushes_type* bushes) {
    int r, i;
    network_type *network = bushes->network;

    deleteVector(bushes->SPcost);
    deleteVector(bushes->flow);
    deleteVector(bushes->nodeFlow);
    deleteVector(bushes->weight);
    deleteVector(bushes->nodeWeight);
    deleteVector(bushes->likelihood);

    for (r = 0; r < network->numZones; r++) {
        for (i = 0; i < network->numNodes; i++) {
            clearArcList(&(bushes->bushForwardStar[r][i]));
            clearArcList(&(bushes->bushReverseStar[r][i]));
        }
    }
    deleteMatrix(bushes->bushForwardStar, network->numZones);
    deleteMatrix(bushes->bushReverseStar, network->numZones);
    deleteVector(bushes->numBushLinks);
    deleteVector(bushes->numBushPaths);
    deleteScalar(bushes);
}

/*
 * bushTopologicalOrder -- Find a topological order using the standard
 * algorithm (finding and marking nodes with no marked predecessors).
 * Arguments are the origin corresponding to the bush, and the network/bush
 * data structures.
 */
void bushTopologicalOrder(int origin, network_type *network,
                          bushes_type *bushes) {
    arcListElt *curArc;
    int i, j, next;
    declareVector(int, indegree, network->numNodes);
    for (i = 0; i < network->numNodes; i++) {
        indegree[i] = bushes->bushReverseStar[origin][i].size;
        bushes->bushOrder[origin][i] = NO_PATH_EXISTS;
    }

    queue_type LIST = createQueue(network->numNodes, network->numNodes);
    next = 0;
    for (i = 0; i < network->numNodes; i++)
        if (indegree[i] == 0)
            enQueue(&LIST, i);
    while (LIST.curelts > 0) {
        i = deQueue(&LIST);
        bushes->bushOrder[origin][next] = i;
        next++;
        for (curArc = bushes->bushForwardStar[origin][i].head;
             curArc != NULL;
             curArc = curArc->next) 
        {
            j = curArc->arc->head;
            indegree[j]--;
            if (indegree[j] == 0) enQueue(&LIST, j);
        }
    }
    if (next < network->numNodes) {
        fatalError("Graph given to bushTopologicalOrder contains a cycle.");
    }

    deleteQueue(&LIST);
    deleteVector(indegree);
}

/* 
 * bushShortestPath -- Given link costs, identify the shortest path
 * using bush links only.  Since the bush is acyclic we can do this
 * very quickly.
 */
void bushShortestPath(network_type *network, bushes_type *bushes, int origin) {
    int curnode, i; /* curnode is topological order, i is real index */
    int h; /* upstream node */
    arcListElt *curArc;

    bushes->SPcost[origin] = 0;
    for (curnode = 1; curnode < network->numNodes; curnode++) {
        i = bushes->bushOrder[origin][curnode];
        bushes->SPcost[i] = INFINITY;
        for (curArc = bushes->bushReverseStar[origin][i].head;
             curArc != NULL;
             curArc = curArc->next)
        {
            h = curArc->arc->tail;
            bushes->SPcost[i] = min(bushes->SPcost[i],
                                    bushes->SPcost[h] + curArc->arc->cost);
        }
    }
}

/*
 * dialFlows -- Use Dial's method to first compute link likelihoods;
 * and then link/node weights; and then link/node flows.
 * These are returned in the flow array of the bushes struct.
 *
 * You need to call bushShortestPath first, or else the likelihoods
 * may not be valid.
 */
void dialFlows(network_type *network, bushes_type *bushes, int origin,
               double theta) {
    int curnode, i, j, ij;
    arcListElt *curArc;

    /* 1. Compute link likelihoods */
    for (ij = 0; ij < network->numArcs; ij++) {
        i = network->arcs[ij].tail;
        j = network->arcs[ij].head;
        bushes->likelihood[ij] = exp(theta * (bushes->SPcost[j]
                                              - bushes->SPcost[i]
                                              - network->arcs[ij].cost));
    }

    /* 2. Compute node/link weights, starting with origin... */
    bushes->nodeWeight[origin] = 1;
    for (curArc = bushes->bushForwardStar[origin][origin].head;
         curArc != NULL;
         curArc = curArc->next)
    {
        ij = ptr2arc(network, curArc->arc);
        bushes->weight[ij] = bushes->likelihood[ij];
    }
    /* ... and now for the other nodes in topological order. */
    for (i = 0; i < network->numNodes; i++) {
        bushes->nodeWeight[i] = 0;
        for (curArc = bushes->bushReverseStar[origin][i].head;
             curArc != NULL;
             curArc = curArc->next)
        {
            ij = ptr2arc(network, curArc->arc);
            bushes->nodeWeight[i] += bushes->weight[ij];
        }
        for (curArc = bushes->bushForwardStar[origin][i].head;
             curArc != NULL;
             curArc = curArc->next)
        {
            ij = ptr2arc(network, curArc->arc);
            bushes->weight[ij] = bushes->nodeWeight[i]*bushes->likelihood[ij];
        }
    }

    /* 3. Now compute node/link flows, in reverse topological order */
    i = bushes->bushOrder[origin][network->numNodes - 1];
    bushes->nodeFlow[i] = (i < network->numZones ?
                           network->demand[origin][i] :
                           0);
    for (curArc = bushes->bushReverseStar[origin][i].head;
         curArc != NULL;
         curArc = curArc->next)
    {
        ij = ptr2arc(network, curArc->arc);
        bushes->flow[ij] = bushes->nodeFlow[i]
                           * (bushes->weight[ij] / bushes->nodeWeight[i]);
    }
    for (curnode = network->numNodes - 2; curnode >= 0; curnode--) {
        bushes->nodeFlow[i] = 0;
        for (curArc = bushes->bushForwardStar[origin][i].head;
             curArc != NULL;
             curArc = curArc->next)
        {
            ij = ptr2arc(network, curArc->arc);
            bushes->nodeFlow[i] += bushes->flow[ij];
        }
        for (curArc = bushes->bushReverseStar[origin][i].head;
             curArc != NULL;
             curArc = curArc->next)
        {
            ij = ptr2arc(network, curArc->arc);
            bushes->flow[ij] = bushes->nodeFlow[i]
                               * (bushes->weight[ij] / bushes->nodeWeight[i]);
        }
    }
}
