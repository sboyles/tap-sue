/*
 * This file contains (1) implementations of general-purpose network algorithms
 * and (2) supporting infrastructure for the network data structure.
 *
 * Regarding (1), this file includes a heap-based implementation of Dijkstra's
 * algorithm, and a network connectivity checker.
 *
 * Regarding (2), this file contains code for displaying network data in
 * human-readable format, and implementations of linked lists for links and
 * paths.
 */

#include "networks.h"

/*
 * Heap-based implementation of Dijkstra's algorithm.  Most of the
 * computations will be done on acyclic networks using a specialized
 * shortest path algorithm, and this will only be called during
 * initialization to generate those bushes by identifying reasonable
 * links.  We only need labels for this, not the paths themselves,
 * so that's all this function does.
 */
void shortestPath(int origin, double *label, network_type *network) {
    int j;
    arcListElt *i;
    int curnode;
   
    /* Initialize heap */
    double tempLabel;
    heap_type *dijkstraHeap = createHeap(network->numNodes, network->numNodes);

    /* Initialize Dijkstra's */
    for (j = 0; j < network->numNodes; j++) {
        dijkstraHeap->valueFn[j] = INFINITY; 
        /* valueFn in the heap stores the cost labels */
    }

    /* Now iterate until the heap is empty */
    insertHeap(dijkstraHeap, origin, 0);
    while (dijkstraHeap->last >= 0) {
        curnode = findMinHeap(dijkstraHeap);
        deleteMinHeap(dijkstraHeap);
        for (i = network->nodes[curnode].forwardStar.head; i != NULL; 
                i = i->next) {
            j = i->arc->head;
            tempLabel = dijkstraHeap->valueFn[curnode] + i->arc->cost;
            if (tempLabel < dijkstraHeap->valueFn[j]) {
                /* Avoid centroid connectors */
                if (j < network->firstThroughNode) {
                    dijkstraHeap->valueFn[j] = tempLabel;
                    continue;
                }

                if (dijkstraHeap->valueFn[j] == INFINITY)
                    insertHeap(dijkstraHeap, j, tempLabel);
                else
                    decreaseKey(dijkstraHeap, j, tempLabel);
            }
        }
    }

    /* Now copy labels to return, and clean up memory */
    memcpy(label, dijkstraHeap->valueFn, sizeof(double) * network->numNodes);
    deleteHeap(dijkstraHeap);
}

/*
finalizeNetwork: After adding the links and nodes to the network struct, this
function generates the forward and reverse star lists. 
*/
void finalizeNetwork(network_type *network) {
    int i, ij, c;

    for (i = 0; i < network->numNodes; i++) {
        initializeArcList(&(network->nodes[i].forwardStar));
        initializeArcList(&(network->nodes[i].reverseStar));
    }
    for (ij = 0; ij < network->numArcs; ij++) {
        insertArcList(&(network->nodes[network->arcs[ij].tail].forwardStar),
                &(network->arcs[ij]), 
                network->nodes[network->arcs[ij].tail].forwardStar.tail);
        insertArcList(&(network->nodes[network->arcs[ij].head].reverseStar),
                &(network->arcs[ij]), 
                network->nodes[network->arcs[ij].head].reverseStar.tail);
        network->arcs[ij].fixedCost = (network->arcs[ij].length
                                          * network->distanceFactor)
                                      + (network->arcs[ij].toll
                                          * network->tollFactor);
        network->arcs[ij].cost = network->arcs[ij].freeFlowTime 
                                 + network->arcs[ij].fixedCost;
        network->arcs[ij].flow = 0;
    }
}

/*
search: Given an initial node (the origin argument), performs a search to
identify all nodes reachable from origin, or from which origin can be reached,
depending on the argument 'd' (FORWARD = nodes reachable from origin; REVERSE =
nodes from which origin can be reached).  Argument 'q' indicates the search
order (FIFO = breadth-first search, LIFO = depth-first search).  Upon
termination, the arrays order and backnode are returned.  backnode indicates
the previous/next node on the path from/to origin (depending on d); if
backnode[i] is the symbolic constant NO_PATH_EXISTS then node i is not
connected.  The order array indicates the order in which nodes are found.
*/
void search(int origin, int* order, int *backnode, network_type *network,
        queueDiscipline q, direction_type d) { int i, j = NO_PATH_EXISTS,
    next; arcListElt *curarc;

    /* Initialize; any node for which backnode[i] remains at NO_PATH_EXISTS is 
     * not connected from/to origin */
    for(i = 0; i < network->numNodes; i++) {
        backnode[i] = NO_PATH_EXISTS;
    }
    backnode[origin] = 0;
    next = 1;
    order[origin] = next;

   /* List of visited nodes is maintained as a queue with discipline q */
    queue_type LIST = createQueue(network->numNodes, network->numNodes);
    enQueue(&LIST, origin);

    /* This code uses a circular queue implementation; queue is empty iff 
     * readptr and writeptr are identical */
    while (LIST.readptr != LIST.writeptr) {
        i = deQueue(&LIST);
        /* Identify the proper list (forward or reverse) ... */
        switch (d) {
            case FORWARD: curarc = network->nodes[i].forwardStar.head; break;
            case REVERSE: curarc = network->nodes[i].reverseStar.head; break;
            default: fatalError("Unknown direction in search."); break;
        }
        /* ...and now iterate through all its elements */
        while (curarc != NULL) {
            switch (d) {
                case FORWARD: j = curarc->arc->head; break;
                case REVERSE: j = curarc->arc->tail; break;
            }
            if (backnode[j] == NO_PATH_EXISTS) { /* Is admissible; arc 
                                                    discovers a new node */
                backnode[j] = i;
                order[j] = ++next;
                if (j >= network->firstThroughNode) {
                    switch (q) {
                        case FIFO: enQueue(&LIST, j); break;
                        case LIFO: frontQueue(&LIST, j); break;
                        case DEQUE:
                            switch (LIST.history[j]) {
                                case NEVER_IN_QUEUE: enQueue(&LIST, j); break;
                                case WAS_IN_QUEUE: frontQueue(&LIST, j); break;
                            }
                            break;
                        default: 
                            fatalError("Unsupported queue type in search."); 
                            break;
                    }
                }
            }
            curarc = curarc->next;
        }
    }
    deleteQueue(&LIST);
    return;
}

/* Update all link costs based on current flows */
void updateLinkCosts(network_type *network) {
    int ij;
    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].cost =
                network->arcs[ij].calculateCost(&network->arcs[ij]);
    }
}

/*
 * generalBPRcost -- Evaluates the BPR function for an arbitrary polynomial.
 */
double generalBPRcost(arc_type *arc) {
   if (arc->flow <= 0)
   // Protect against negative flow values and 0^0 errors
       return arc->freeFlowTime + arc->fixedCost;

   return arc->fixedCost + arc->freeFlowTime *
       (1 + arc->alpha * pow(arc->flow / arc->capacity, arc->beta));
}

/* linearBPRcost -- Faster implementation for linear BPR functions. */
double linearBPRcost(arc_type *arc) {
   return arc->fixedCost + arc->freeFlowTime *
       (1 + arc->alpha * arc->flow / arc->capacity);
}

/* quarticBPRcost -- Faster implementation for 4th-power BPR functions
 */
double quarticBPRcost(arc_type *arc) {
   double y = arc->flow / arc->capacity;
   y *= y;
   y *= y;
   return arc->fixedCost + arc->freeFlowTime * (1 + arc->alpha * y);
}

/*
forwardStarOrder is a comparison function; a pointer to this function can be
passed to qsort or other sorting routines.  In forward star order, a link
precedes another if its tail node has a lower index.  This function implements
a tiebreaking rule based on the head node.
*/
int forwardStarOrder(const void *arc1, const void *arc2) {
    arc_type first = *(arc_type *)arc1; 
    arc_type second = *(arc_type *)arc2;
    if (first.tail <  second.tail) {
        return -1;
    } else if (first.tail > second.tail) {
        return 1;
    } else if (first.head < second.head) {
        return -1;
    } else if (first.head > second.head) {
        return 1;
    } else { 
        return 0;
    }
}

/*
ptr2arc converts a pointer to an arc to the index number for that arc; to do
this, the network struct needs to be passed aint with the arc pointer.
*/
int ptr2arc(network_type *network, arc_type *arcptr) {
    return (int) (arcptr - network->arcs);
}

/*
displayNetwork prints network data in human-readable format.  minVerbosity is
used to control whether anything needs to be printed.
*/

void displayNetwork(int minVerbosity, network_type *network) {
    int i;
    displayMessage(minVerbosity, "Network has %d nodes and %d arcs\n", 
            network->numNodes, network->numArcs);
    displayMessage(minVerbosity, "Arc data: ID, tail, head, flow, cost, der "
                                 "(skipping artificial arcs)\n");
    for (i = 0; i < network->numArcs; i++) {
       if (network->arcs[i].capacity == ARTIFICIAL) continue; 
       displayMessage(minVerbosity, "%ld (%ld,%ld) %f %f %f\n", i, 
               network->arcs[i].tail + 1, network->arcs[i].head + 1, 
               network->arcs[i].flow, network->arcs[i].cost, 
               network->arcs[i].der);
    }
}

/*
deleteNetwork deallocates any memory assigned to a network struct.
*/
void deleteNetwork(network_type *network) {
   int i;
   for (i = 0; i < network->numNodes; i++) {
      clearArcList(&(network->nodes[i].forwardStar));
      clearArcList(&(network->nodes[i].reverseStar));
   }
   for (i = 0; i < network->numArcs; i++) {
       deleteVector(network->arcs[i].classFlow);
       deleteVector(network->arcs[i].classCost);
       deleteVector(network->arcs[i].classToll);
   }

   deleteMatrix(network->demand, network->batchSize);
   deleteVector(network->nodes);
   deleteVector(network->arcs);
   deleteVector(network->tollFactor);
   deleteVector(network->distanceFactor);
   deleteScalar(network);
}

/*
The functions below implement doubly-linked lists of arcs.
You probably don't need to poke around here too much.
*/
arcList *createArcList() {
    declareScalar(arcList, newdll);
    initializeArcList(newdll);
    return newdll;
}

void initializeArcList(arcList *list) {
    list->head = NULL;
    list->tail = NULL;
    list->size = 0;
}

arcListElt *insertArcList(arcList *list, arc_type *value, arcListElt *after) {
    declareScalar(arcListElt, newNode);
    newNode->arc = value;
    if (after != NULL) {
        newNode->prev = after;
        newNode->next = after->next;
        if (list->tail != after) 
            newNode->next->prev = newNode; 
        else 
            list->tail = newNode;
        after->next = newNode;
    } else {
        newNode->prev = NULL;
        newNode->next = list->head;
        if (list->tail != after) 
            newNode->next->prev = newNode; 
        else 
            list->tail = newNode;
        list->head = newNode;
    }
    list->size++;
    return newNode;
}

void clearArcList(arcList *list) {
    while (list->head != NULL)
        deleteArcListElt(list, list->tail);
}

void deleteArcList(arcList *list) {
    clearArcList(list);
    deleteScalar(list);
}

void deleteArcListElt(arcList *list, arcListElt *elt) {
    if (list->tail != elt) {
        if (list->head != elt) 
            elt->prev->next = elt->next; 
        else 
            list->head = elt->next;
        elt->next->prev = elt->prev;
    } else {
        list->tail = elt->prev;
        if (list->head != elt) 
            elt->prev->next = elt->next; 
        else 
            list->head = elt->next;
    }
    list->size--;
    deleteScalar(elt);
}

void displayArcList(arcList *list) {
    arcListElt *curnode = list->head;
    printf("Start of the list: %p\n", (void *)list->head);
    while (curnode != NULL) {
        printf("%p (%d,%d) %p %p\n", (void *)curnode, curnode->arc->tail, 
                curnode->arc->head, (void *)curnode->prev, 
                (void *)curnode->next);
        curnode = (*curnode).next;
    }
    printf("End of the list: %p\n", (void *)list->tail);
}
