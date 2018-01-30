/* options.h: heder file for structure Options
 */

#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#include "util.h"
#include "algor.h"

/* ===========================================================================
 *                              structure Options definition
 * ===========================================================================
 */
typedef struct {
    uint32_t seed;
    int32_t algor;                  /* which algorithm to use */
    
    uint64_t crossbred_mb_size;     /* number of solutions to keep in the mailbox */
    uint64_t crossbred_kvar_num;    /* number of variables to keep */
    uint64_t crossbred_mac_degree;  /* degree of Macaulay matrix to create */
    uint64_t crossbred_kf_num;      /* number of variables to fix by kernel */
    uint64_t crossbred_mac_keq_num; /* number of equations to keep for
                                     * Macaulay matrix reduction
                                     */
    uint64_t crossbred_sub_keq_num; /* number of equations to keep in the sub-system */
    uint64_t crossbred_sub_fnum;    /* number of variables to fix after a
                                     * sub-system is obtained
                                     */
    double crossbred_mq_ratio;      /* the ratio of monomials present in the initial
                                     * MQ system
                                     */
    uint64_t thread_num;            /* number of CPU threads to use */
    char* cha_file;
    char* mq_fix_file;              /* varibles to fix before applying anything */

    int crossbred_dev_id;           /* specify which GPU device to use */
    bool crossbred_mac_stats;       /* if the program should collects stats
                                     * about Macaulay matrix
                                     */
    bool verbose;
    bool crossbred_rmac_cpu;        /* use cpu to perform Gaussian elimination
                                     * on reduced Macaulay matrix instead of gpu
                                     */
    bool crossbred_init_gj;         /* perform Gaussian-Jordan elimination on the
                                     *  input MQ system. Default to true.
                                     */
    uint32_t crossbred_mode;        /* operation mode, see macros defined below */
    char* mq_ext_file;              /* file for storing the extracted sub-systems */
} Options;


#define MQ_MODE_DEFAULT 0           /* extract sub-systems and bruteforce */
#define MQ_MODE_EXT     1           /* only extract a sub-system */
#define MQ_MODE_BF      2           /* read an extracted sub-system and solve it */

/* ===========================================================================
 *                              function prototypes
 * ===========================================================================
 */

/* function: options_init
 * usage: init the Options structure
 * arguments: an pointer to structure Options
 * return: void
 */
void
options_init(Options* opts);

/* function: options_free
 * usage: free the structure Options
 * arguments an pointer to structure Options
 * return: void
 */
void
options_free(Options* opts);

/* function: options_parse
 * usage: parse argv and store the options into structure Options
 *      Make sure options_init has been called.
 * arguments:
 *      1) opts: an pointer to the Options structure
 *      2) argc
 *      3) argv
 * return: void
 */
void
options_parse(Options* opts, int argc, char** argv);

#endif /* __OPTIONS_H__ */
