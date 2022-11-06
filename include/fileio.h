/*
 * fileio.h -- This is the header for file reading and writing.  String
 * processing routines also go here.
 */

#ifndef FILEIO_H
#define FILEIO_H

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bush.h"
#include "networks.h"
#include "utils.h"

#define STRING_SIZE 9999

#ifdef DEBUG_MODE
char debugFileName[STRING_SIZE];
FILE *debugFile;
#endif

enum { // Return codes for metadata parsing
    SUCCESS,
    BLANK_LINE,
    COMMENT
};

///////////////////////////
// Reading network files //
///////////////////////////

void readTntpNetwork(network_type *network, char *linkFileName,
                    char *tripFileName);

///////////////////////
// String processing //
///////////////////////

void blankInputString(char *string, int length);
int parseMetadata(char* inputLine, char* metadataTag, char* metadataValue);
int parseLine(char* inputLine, char* outputLine);
#endif
