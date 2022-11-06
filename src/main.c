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
    network_type *network = newScalar(network_type);
#ifdef DEBUG_MODE
    debugFile = openFile("full_log.txt", "w");
#endif

    if (argc != 5)
        fatalError("Must specify exactly four parameters "
                   "(network file, trips file, theta, lambda).\n");
    readTntpNetwork(network, argv[1], argv[2]);
    SUE_MSA(network, atof(argv[3]), atof(argv[4]));
    deleteNetwork(network);

#ifdef DEBUG_MODE
    fclose(debugFile);
#endif

    return EXIT_SUCCESS;
}
