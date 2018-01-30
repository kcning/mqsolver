/* main.c: main program of the MQ challenge solver
 */

#include "mqsolver.h"


int
main(int argc, char* argv[]) {
    MQSolver mqsolver;
    mqs_init(&mqsolver, argc, argv);
    mqs_solve(&mqsolver);
    mqs_free(&mqsolver);
    return 0;
}
