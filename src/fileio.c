#include "fileio.h"

///////////////////////////
// Reading network files //
///////////////////////////

void readTntpNetwork(network_type *network, char *linkFileName, 
                     char *tripFileName) {
    int i, j, r = 0;
    int check;
    int numParams, status;
    double demand, totalDemandCheck = 0;
    double defaultTollFactor, defaultDistanceFactor;

    char fullLine[STRING_SIZE], trimmedLine[STRING_SIZE], *token;
    char metadataTag[STRING_SIZE], metadataValue[STRING_SIZE];

    FILE *linkFile = openFile(linkFileName, "r");
    FILE *tripFile;

    network->numZones = IS_MISSING;
    network->numArcs = IS_MISSING;
    network->numNodes = IS_MISSING;
    network->firstThroughNode = IS_MISSING;
    defaultTollFactor = IS_MISSING;
    defaultDistanceFactor = IS_MISSING;

    /* Read link file metadata */
    bool endofMetadata = FALSE;
    do {
        if (fgets(fullLine, STRING_SIZE, linkFile) == NULL) 
            fatalError("Link file %s ended (or other I/O error) before "
                        "metadata complete.", linkFileName);
        status = parseMetadata(fullLine, metadataTag, metadataValue);
        if (status == BLANK_LINE || status == COMMENT) continue;
        if         (strcmp(metadataTag, "NUMBER OF ZONES") == 0) {
            network->numZones = atoi(metadataValue);
        } else if (strcmp(metadataTag, "NUMBER OF LINKS") == 0) {
            network->numArcs = atoi(metadataValue);
        } else if (strcmp(metadataTag, "NUMBER OF NODES") == 0) {
            network->numNodes = atoi(metadataValue);
        } else if (strcmp(metadataTag, "FIRST THRU NODE") == 0) {
            network->firstThroughNode = atoi(metadataValue) - 1;
        } else if (strcmp(metadataTag, "DISTANCE FACTOR") == 0) {
            defaultDistanceFactor = atof(metadataValue);
        } else if (strcmp(metadataTag, "TOLL FACTOR") == 0) {
            defaultTollFactor = atof(metadataValue);
        } else if (strcmp(metadataTag, "END OF METADATA") == 0) {
            endofMetadata = TRUE;
        } else {
            warning(MEDIUM_NOTIFICATIONS, "Ignoring unknown metadata tag %s "
                    "in link file %s", metadataTag, linkFileName);
        }
    } while (endofMetadata == FALSE);

    /* Check input for completeness and correctness */
    if (network->numZones == IS_MISSING) 
        fatalError("Link file %s does not contain number of zones.", 
                linkFileName);
    if (network->numNodes == IS_MISSING) 
        fatalError("Link file %s does not contain number of nodes.", 
                linkFileName);
    if (network->numArcs == IS_MISSING) 
        fatalError("Link file %s does not contain number of links.", 
                linkFileName);
    if (network->firstThroughNode == IS_MISSING) {
        warning(LOW_NOTIFICATIONS, "Link file %s does not contain first "
                "through node, setting to 1 as default.\n", linkFileName);
        network->firstThroughNode = 0;
    }
    if (defaultDistanceFactor == IS_MISSING) {
        defaultDistanceFactor = 0;
    }
    if (defaultTollFactor == IS_MISSING) {
        defaultTollFactor = 0;
    }
    if (network->numZones < 1) 
        fatalError("Link file %s does not contain a positive number of nodes.",
                linkFileName);
    if (network->numArcs < 1) 
        fatalError("Link file %s does not contain a positive number of links.",
                linkFileName);
    if (network->numNodes < 1) 
        fatalError("Link file %s does not contain a positive number of nodes.",
                linkFileName);

    displayMessage(MEDIUM_NOTIFICATIONS, "Nodes, arcs, zones, thrunode: "
            "%ld %ld %ld %ld\n", network->numNodes, network->numArcs,
            network->numZones, network->firstThroughNode);
    displayMessage(MEDIUM_NOTIFICATIONS, "Distance factor, toll factor: "
            "%lf %lf\n", defaultDistanceFactor, defaultTollFactor);

    network->nodes = newVector(network->numNodes, node_type);
    network->arcs = newVector(network->numArcs, arc_type);
    network->demand = newMatrix(network->numZones, network->numZones,double);
    for (i = 0; i < network->numZones; i++) {
        for (j = 0; j < network->numZones; j++) {
            network->demand[i][j] = 0;
        }
    }

    /* Read link data */
    for (i = 0; i < network->numArcs; i++) {
        if (fgets(fullLine, STRING_SIZE, linkFile) == NULL)
            fatalError("Link file %s ended (or other I/O error) before link "
                    "data complete.", linkFileName);
        status = parseLine(fullLine, trimmedLine);
        if (status == BLANK_LINE || status == COMMENT) {
            i--;
            continue;
        }
        numParams=sscanf(trimmedLine,"%d %d %lf %lf %lf %lf %lf %lf %lf %d",
            &network->arcs[i].tail,
            &network->arcs[i].head,
            &network->arcs[i].capacity,
            &network->arcs[i].length,
            &network->arcs[i].freeFlowTime,
            &network->arcs[i].alpha,
            &network->arcs[i].beta,
            &network->arcs[i].speedLimit,
            &network->arcs[i].toll,
            &network->arcs[i].linkType);
        if (numParams != 10) 
            fatalError("Link file %s has an error in this line:\n\"%s\"",
                    linkFileName, fullLine);
        if (network->arcs[i].tail < 1 
                || network->arcs[i].tail > network->numNodes) 
            fatalError("Arc tail %d out of range in network file %s.", 
                    i, linkFileName);
        if (network->arcs[i].head < 1 
                || network->arcs[i].head > network->numNodes) 
            fatalError("Arc head %d out of range in network file %s.", 
                    i, linkFileName);
        if (network->arcs[i].length < 0) 
            warning(FULL_NOTIFICATIONS, 
                    "Arc length %d negative in network file %s.\n%s", i,
                    linkFileName, fullLine);
        if (network->arcs[i].freeFlowTime < 0) 
            fatalError("Arc free flow time %d negative in network file "
                    "%s.\n%s", i, linkFileName, fullLine);
        if (network->arcs[i].alpha < 0) fatalError("Alpha %d negative in "
                "network file %s.\n%s", i, linkFileName, fullLine);
        if (network->arcs[i].beta < 0) fatalError("Beta %d negative in "
                "network file %s.\n%s", i, linkFileName, fullLine);
        if (network->arcs[i].speedLimit < 0) warning(FULL_NOTIFICATIONS, 
                "Speed limit %d negative in network file %s.\n%s", i, 
                linkFileName, fullLine);
        if (network->arcs[i].toll < 0) warning(FULL_NOTIFICATIONS, "Toll %d "
                "negative in network file %s.\n%s", i, linkFileName, fullLine);
        if (network->arcs[i].capacity <= 0) fatalError("Capacity %d "
                "nonpositive in network file %s.\n%s", i, linkFileName, 
                fullLine);
        network->arcs[i].tail--;
        network->arcs[i].head--;
        network->arcs[i].flow = 0;
        network->arcs[i].cost = network->arcs[i].freeFlowTime;
        if (network->arcs[i].beta == 1) {
           network->arcs[i].calculateCost = &linearBPRcost;
        } else if (network->arcs[i].beta == 4) {
           network->arcs[i].calculateCost = &quarticBPRcost;
        } else {
           network->arcs[i].calculateCost = &generalBPRcost;
        }           
    }
    fclose(linkFile);

    tripFile = openFile(tripFileName, "r");
    /* Verify trip table metadata */
    endofMetadata = FALSE;
    network->totalODFlow = IS_MISSING;
    network->tollFactor = defaultTollFactor;
    network->distanceFactor = defaultDistanceFactor;
    do {
        if (fgets(fullLine, STRING_SIZE, tripFile) == NULL) 
            fatalError("Trip file %s ended (or other I/O error) before "
                    "metadata complete.", tripFileName);
        status = parseMetadata(fullLine, metadataTag, metadataValue);
        if (status == BLANK_LINE || status == COMMENT) continue;
        if         (strcmp(metadataTag, "NUMBER OF ZONES") == 0) {
            check = atoi(metadataValue);
            if (check != network->numZones) fatalError("Number of zones in"
                    "trip and link files do not match.");
        } else if (strcmp(metadataTag, "TOTAL OD FLOW") == 0) {
            network->totalODFlow = atof(metadataValue);
        } else if (strcmp(metadataTag, "DISTANCE FACTOR") == 0) {
            network->distanceFactor = atof(metadataValue);
        } else if (strcmp(metadataTag, "TOLL FACTOR") == 0) {
            network->tollFactor = atof(metadataValue);
        } else if (strcmp(metadataTag, "END OF METADATA") == 0) {
            endofMetadata = TRUE;
        } else {
            warning(MEDIUM_NOTIFICATIONS, "Ignoring unknown metadata tag "
                    "%s in trips file %s", metadataTag, tripFileName);
        }
    } while (endofMetadata == FALSE);

    /* Now read trip table */
    while (!feof(tripFile)) {
        if (fgets(fullLine, STRING_SIZE, tripFile) == NULL) break;
        status = parseLine(fullLine, trimmedLine);
        if (status == BLANK_LINE || status == COMMENT) continue;
        if (strstr(trimmedLine, "Origin") != NULL) {
            // i indexes current origin
            sscanf(strstr(trimmedLine, "Origin")+6,"%d", &i);  
            if (i <= 0 || i > network->numZones) fatalError("Origin %d is"
                    "out of range in trips file %s", i, tripFileName);
            i--;
            continue;
        }
        token = strtok(trimmedLine , ";");
        while (token != NULL && strlen(token) > 1) {
            numParams = sscanf(token, "%d : %lf", &j, &demand);
            if (numParams < 2) break;
            if (j <= 0 || j > network->numZones) fatalError("Destination "
                    "%d is out of range in trips file %s\n%s\n%s", j, 
                    tripFileName, fullLine, token);
            j--;
            network->demand[r][j] = demand;
            if (demand < 0) fatalError("Negative demand from origin %d to "
                    "destination %d", i, j);
            totalDemandCheck += network->demand[i][j];
            token = strtok(NULL, ";");
        }
        blankInputString(trimmedLine, STRING_SIZE);
    }

    fclose(tripFile);
    finalizeNetwork(network);
    displayMessage(FULL_NOTIFICATIONS, "Forward and reverse star lists "
            "generated.\n");
}

///////////////////////
// String processing //
///////////////////////

void blankInputString(char *string, int length) {
    int i;
    for (i = 0; i < length; i++) string[i] = '\0';
}

int parseMetadata(char* inputLine, char* metadataTag, char* metadataValue) {
    /* metadataTag and metadataValue both need to be of at least STRING_SIZE */
    int i = 0, j = 0;
    inputLine[STRING_SIZE-1] = '\0';
    while (inputLine[i] != '\0' && inputLine[i] != '\n' && inputLine[i] != '\r' 
            && inputLine[i] != '<' && inputLine[i] != '~') i++;
    if (inputLine[i] == '\0' || inputLine[i] == '\n' || inputLine[i] == '\r') 
        return BLANK_LINE;
    if (inputLine[i] == '~') return COMMENT;
    i++;
    while (inputLine[i] != '\0' && inputLine[i] != '>') {
        metadataTag[j++] = toupper(inputLine[i++]);
    }
    metadataTag[j] = '\0';
    if (inputLine[i] == '\0')
        fatalError("Metadata tag not closed: ", metadataTag);
    i++;
    while (inputLine[i] != '\0' 
            && (inputLine[i] == ' ' || inputLine[i] == '\t')) 
        i++;
    j = 0;
    while (inputLine[i] != '\0' && inputLine[i] != '\n' && inputLine[i] != '~')
        metadataValue[j++] = inputLine[i++];
    metadataValue[j] = '\0';
    return SUCCESS;
}

// Checks for comments and blank lines, and removes leading spaces
int parseLine(char* inputLine, char* outputLine) {
    int i = 0, j = 0;
    while (inputLine[i] != '\0' && (inputLine[i] == ' ' 
                || inputLine[i] == '\t')) i++;
    if (inputLine[i] == '~') return COMMENT;
    if (inputLine[i] == '\0' || inputLine[i] == '\n' 
            || inputLine[i] == '\r') return BLANK_LINE;
    while (inputLine[i] != '\0' && i < STRING_SIZE - 1) {
        outputLine[j++] = inputLine[i++];
    }
    outputLine[j] = '\0';
    return SUCCESS;
}
