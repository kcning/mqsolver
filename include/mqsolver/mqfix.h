/* mqfix.h: heder file for structure MQfix
 */

#ifndef __MQFIX_H__
#define __MQFIX_H__

#include <stdbool.h>
#include <stdint.h>

#define MQFIX_SEQUENTIAL        0x0UL
#define MQFIX_FILE              0x1UL

typedef struct {
    uint64_t fnum;
    uint64_t lnum;
    uint64_t type;

    /* for fixing variables in a sequential manner */
    uint64_t start;  /* included */
    uint64_t end;  /* not included */

    /* for fixing variables arbitrarily */
    uint64_t val[];  /* each line is stored as a uint64_t therefore the maximal
                      * variables that can be fixed during preprocessing is 64.
                      * x_{n-f+1} is stored as the last bit and so on.
                      */
} MQFix;

/* function: mqfix_create
 * usage: create an MQFix struct with a file that specifies how
 *      to fix variables in the original MQ system
 * arguments:
 *      1) file: path to the file
 * return: address to the newly created structure MQFix. On error, NULL is
 *      returned
 */
MQFix*
mqfix_create(const char* file);

/* function: mqfix_create_seq
 * usage: create an MQFix struct by specifying the starting point and the range
 * arguments:
 *      1) fnum: number of variables to fix during preprocessing
 *      2) start: the first preprocessing configuration
 *      3) range: the range within which preprocessing configs should be used
 * return: address to the newly created structure MQFix. On error, NULL is
 *      returned
 */
MQFix*
mqfix_create_seq(const uint64_t fnum, const uint64_t start,
                 const uint64_t range);

/* function: mqfix_free
 * usage: free the MQFix structure and its fields
 * arguments:
 *      1) a pointer to structure MQFix
 * return: void
 */
void
mqfix_free(MQFix* mqf);

/* function: mqfix_line
 * usage: return one line, should only be called when a configuration file
 *      is used to create MQFix
 * arguments:
 *      1) a pointer to structure MQFix
 *      2) lidx: index of the line
 * return: address to the line
 */
uint64_t*
mqfix_line(MQFix* mqf, const uint64_t lidx);

/* function: mqfix_lvar
 * usage: given index to a line and the index of the target varible, return the
 *      value for the variable
 * arguments:
 *      1) a pointer to structure MQFix
 *      2) lidx: index to the line
 *      3) idx: index of the variable. For x_{n-f+1}, pass 0, and so on
 * return: the value for that variable
 */
bool
mqfix_lvar(MQFix* mqf, const uint64_t lidx, const uint64_t idx);

/* function: mqfix_ltob
 * usage: given index to a line, convert it to an array of bool
 * arguments:
 *      1) a pointer to structure MQFix
 *      2) lidx: index to the line
 *      3) dst: container for the bool array, must be at least as large as the
 *              number of fnum
 * return: void
 */
void
mqfix_ltob(MQFix* mqf, const uint64_t lidx, bool* const dst);

#endif /* __MQFIX_H__ */
