/*
 * Main function for solving stochastic user assignment using the
 * method of successive averages.
 */

#include "main.h"

int main(int argc, char* argv[]) {
    /* verbosity is a global variable controlling how much output to produce,
     * see utils.h for possible values*/
    verbosity = FULL_DEBUG;
    //verbosity = FULL_NOTIFICATIONS;
    network_type *network = createScalar(network_type);
#ifdef DEBUG_MODE
    debugFile = openFile("full_log.txt", "w");
#endif

    if (argc != 5)
        fatalError("Must specify exactly four parameters "
                   "(network file, trips file, theta, lambda).\n");
    network = readOBANetwork(network, argc[1], argc[2]);
    SUE_MSA(network, atof(argc[3]), atof(argc[4]));
    deleteNetwork(network);

#ifdef DEBUG_MODE
    fclose(debugFile);
#endif

    return EXIT_SUCCESS;
}
