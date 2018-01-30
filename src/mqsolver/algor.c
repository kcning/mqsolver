/* algor.c: implementation for algor.h
 */

#include "algor.h"

/* function: algor2code
 * usage: parse algorithm string and return corresponding integer code. Return
 *      -1 for invalid argument.
 * arguments: arg, argorithm string passed by the user. Valid strings are:
 *              fast_ex: fast exhaustive search
 *              crossbred: Joux and Vitse's algorithm
 * return: the corresponding integer code
 */

/* valid string identifiers for algorithms */
static const char* mq_alg_fast_ex = "fast_ex";
static const char* mq_alg_crossbred = "crossbred";

/* comparing the algorithm argument */
static inline bool
cmp_alg(const char* arg, const char* alg_str) {
    return (0 == strncmp(arg, alg_str, strnlen(arg, MQ_ALG_MAX_LEN)));
}

int
algor2code(const char* str) {
    if(cmp_alg(str, mq_alg_fast_ex)) {
        return MQ_ALG_FAST_EX;
    } else if(cmp_alg(str, mq_alg_crossbred)) {
        return MQ_ALG_CROSSBRED;
    } else {
        return MQ_ALG_INVALID;
    }
}
