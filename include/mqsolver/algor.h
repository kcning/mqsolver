/* algor.h: define shared code for algorithms for solving MQ */

#ifndef __ALGOR_H__
#define __ALGOR_H__

#include <stdbool.h>
#include <string.h>

#define MQ_ALG_MAX_LEN                  30  /* max len of algorithm argument */

/* code for algorithms */
#define MQ_ALG_INVALID                  (-1)
#define MQ_ALG_FAST_EX                  0
#define MQ_ALG_CROSSBRED                1

/* function: algor2code
 * usage: parse algorithm string and return corresponding integer code. Return
 *      -1 for invalid argument.
 * arguments: arg, argorithm string passed by the user. Valid strings are:
 *              fast_ex: fast exhaustive search
 *              crossbred: Joux and Vitse's algorithm
 * return: the corresponding integer code
 */
int
algor2code(const char* str);

#endif  /* __ALGOR_H__ */
