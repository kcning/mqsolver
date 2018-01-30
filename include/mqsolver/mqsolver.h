/* mqsolver.h: header file for structure MQSolver */

#ifndef __MQSOLVER_H__
#define __MQSOLVER_H__

#include <stdio.h>
#include <string.h>
#include <limits.h>


#include "util.h"
#include "options.h"
#include "algor.h"
#include "mq_math.h"
#include "macaulay.h"
#include "fix.h"
#include "cuda_util.h"
#include "graycode.h"
#include "debug.h"
#include "threadpool.h"
#include "rmac.h"

/* ===========================================================================
 *                              structure MQSolver definition
 * ===========================================================================
 */

typedef struct {
    Options* opts;

    uint64_t var_num;  // number of variables/unknowns of the original system
    uint64_t eq_num;  // number of equations of the original system
    int seed;  // seed in the challenge file
    
    /* number of extended variables, including the original variables,
     * quadratic products of them, and the constant term, e.g. for x1, x2,
     * the set of extended variables contains: x1^2, x1x2, x2^2, x1, x2, 1
     */
    uint64_t xvar_num;

    /* equations */
    bool** orig_sys;
} MQSolver;

/* ===========================================================================
 *                              function prototypes
 * ===========================================================================
 */

/* function: mqsolver_init
 * usage: init the MQSolver structure
 *      1) parse options
 *      2) read system info
 *      3) init random number generator
 *      4) read the challenge file passed via Options
 * arguments:
 *      1) mqs: a pointer to structure MQSolver
 *      2) argc
 *      3) argv
 * return: void
 */
void
mqs_init(MQSolver* mqs, int argc, char* argv[]);

/* function: mqs_free
 * usage: free the MQSolver structure and its fields
 * arguments: a pointer to structure MQSolver
 * return: void
 */
void
mqs_free(MQSolver* mqs);

/* function: mqs_print_orig_sys
 * usage: print the original MQ system to terminal
 * arguments: a pointer to structure MQSolver
 * return: void
 */
void
mqs_print_orig_sys(MQSolver* mqs);

/* function: mqs_solve
 * usage: try to solve the MQ system. The algorithm to be used depends on the
 *       options passed via command line. Once a solution is found, the
 *       function prints it to the terminal and returns.
 * arguments: a pointer to structure MQSolver
 * return: 0 is a solution is found, otherwise 1.
 */
int
mqs_solve(MQSolver* mqs);

#endif /* __MQSOLVER_H__ */
