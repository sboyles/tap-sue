/*
 * networks.h -- header file for general network management.  Defines the data
 * structures for storing networks and their components (arcs, nodes, OD
 * pairs).  Also contains implementations of standard network algorithms
 * (shortest path, connectivity) and routines for managing the network
 * data structures.
 */

#ifndef NETWORKS_H
#define NETWORKS_H

#include <math.h>
#include <string.h>
#include "datastructures.h"

#define NO_PATH_EXISTS -1
#define ARTIFICIAL 99999 /* Value used for costs, etc. on artificial links
                            generated to ensure strong connectivity */

typedef enum {
    FORWARD,
    REVERSE
} direction_type;

/*
 * arc_type -- struct for storing all information for arcs.  Most of the
 * members are self-explanatory (tail, head, flow, cost, etc.) and store data
 * from TNTP files.  
 *
 * Note the calculateCost function pointer which can be tailored to a specific
 * function type for greater efficiency.
 */
typedef struct arc_type {
    int    tail;
    int    head;
    double  flow;
    double  cost;

    /* Main link data */
    double  freeFlowTime;
    double  capacity;
    double  length;
    double  toll;

    /* BPR values */
    double  alpha; 
    double  beta;

    /* Other data provided in TNTP format */
    double  speedLimit;
    int     linkType;

    double  fixedCost; /* Reflects toll and distance */
    double  (*calculateCost)(struct arc_type *arc);
} arc_type;


/* Data structures for linked lists of arcs (used for forward/reverse stars) */
#define arcListElt struct AL
arcListElt {
    arc_type    *arc;
    arcListElt  *prev;
    arcListElt  *next;
};

typedef struct {
    arcListElt  *head;
    arcListElt  *tail;
    int         size;
} arcList;

typedef struct {
    arcList *arcs;
    double  cost;
    double  der;
} path_type;

/* node_type -- data structure for nodes.  Only contains lists of arcs entering
 * and leaving the node */
typedef struct {
    arcList forwardStar;
    arcList reverseStar;
} node_type;


/* network_type -- data structure for the entire network, including arrays of
 * nodes, arcs, and OD pairs, and network size information.  The beckmann and
 * beckmannLB members are used in certain gap calculations.
 */
typedef struct {
    node_type* nodes;
    arc_type*  arcs;
    double**   demand;
    int numNodes;
    int numArcs;
    int numZones; 
    int firstThroughNode;
    double totalODFlow;
    double tollFactor;
    double distanceFactor;
} network_type;

void shortestPath(int origin, double *label, network_type *network);
void finalizeNetwork(network_type *network);
void search(int origin, int* order, int *backnode, network_type *network,
            queueDiscipline q, direction_type d);

void updateLinkCosts(network_type *network);
double generalBPRcost(arc_type *arc);
double linearBPRcost(arc_type *arc);
double quarticBPRcost(arc_type *arc);

int forwardStarOrder(const void *arc1, const void *arc2);
int ptr2arc(network_type *network, arc_type *arcptr);

void deleteNetwork(network_type *network);
void displayNetwork(int minVerbosity, network_type *network);

/////////////////////////////////
// Custom linked lists for TAP //
/////////////////////////////////

arcList *createArcList();
void initializeArcList(arcList *list);
arcListElt *insertArcList(arcList *list, arc_type *value, arcListElt *after);
void clearArcList(arcList *list);
void deleteArcList(arcList *list);
void deleteArcListElt(arcList *list, arcListElt *elt);
void displayArcList(arcList *list);
#endif

