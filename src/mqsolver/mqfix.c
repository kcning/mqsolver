/* mqfix.c: implementation of structure MQFix
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>


#include "util.h"
#include "mqfix.h"

/* function: mqfix_create
 * usage: create an MQFix struct with a file that specifies how
 *      to fix variables in the original MQ system
 * arguments:
 *      1) file: path to the file
 * return: address to the newly created structure MQFix. On error, NULL is
 *      returned
 */
MQFix*
mqfix_create(const char* file) {
    FILE* fp = NULL;
    const size_t buf_size = 0x1 << 20; // 1MB per line
    char* buf = SMALLOC(char, buf_size);

    if(NULL == (fp = fopen(file, "r"))) {
        PRINTF_ERR_STAMP("[!] error opening file: %s\n", file);
        return NULL;
    }
    
    if(NULL == fgets(buf, buf_size, fp)) { // read the first line
        PRINTF_ERR_STAMP("[!] error reading from file: %s\n", file);
        fclose(fp);
        return NULL;
    }

    char* ptr = strtok(buf, " \n");
    if(!ptr) {
        PRINTF_ERR_STAMP("[!] error parsing file: %s\n", file);
        fclose(fp);
        return NULL;
    }
    uint64_t fnum = strtoul(ptr, NULL, 0);

    ptr = strtok(NULL, " \n");
    if(!ptr) {
        PRINTF_ERR_STAMP("[!] error parsing file: %s\n", file);
        fclose(fp);
        return NULL;
    }
    uint64_t lnum = strtoul(ptr, NULL, 0);

    ptr = strtok(NULL, " \n");
    if(ptr) { // sequential mode
        uint64_t start = strtoul(ptr, NULL, 0);
        free(buf);
        return mqfix_create_seq(fnum, start, lnum);
    }

    MQFix* mqf = (MQFix* ) malloc(sizeof(MQFix) + sizeof(uint64_t) * lnum);
    if(!mqf) {
        PRINTF_ERR_STAMP("[!] insufficient memory\n");
        fclose(fp);
        return NULL;
    }
    memset(mqf, 0x0, sizeof(MQFix) + sizeof(uint64_t) * lnum);

    mqf->fnum = fnum;
    mqf->lnum = lnum;
    mqf->start = 0;
    mqf->end = 0;
    mqf->type = MQFIX_FILE;

    for(uint64_t i = 0; i < lnum; ++i) {
        if(NULL == fgets(buf, buf_size, fp)) {
            PRINTF_ERR_STAMP("[!] error parsing file: %s\n", file);
            fclose(fp);
            free(mqf);
            return NULL;
        }

        ptr = strtok(buf, " \n");
        uint64_t vidx = 0;
        while(NULL != ptr) {
            *mqfix_line(mqf, i) |= strtoul(ptr, NULL, 0) << vidx++;
            ptr = strtok(NULL, " \n");
        }

        if(vidx != fnum) {
            PRINTF_ERR_STAMP("[!] inconsistent configuration: %s\n", file);
            fclose(fp);
            free(mqf);
            return NULL;
        }
    }

    fclose(fp);
    free(buf);
    return mqf;
}

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
                 const uint64_t range) {
    MQFix* mqf = (MQFix* ) malloc(sizeof(MQFix));
    if(!mqf) {
        PRINTF_ERR_STAMP("[!] insufficient memory\n");
        return NULL;
    }

    mqf->fnum = fnum;
    mqf->lnum = range;
    mqf->start = start;
    mqf->end = start + range;
    mqf->type = MQFIX_SEQUENTIAL;

    return mqf;
}

/* function: mqfix_free
 * usage: free the MQFix structure and its fields
 * arguments:
 *      1) a pointer to structure MQFix
 * return: void
 */
inline void
mqfix_free(MQFix* mqf) {
    free(mqf);
}

/* function: mqfix_line
 * usage: return one line, should only be called when a configuration file
 *      is used to create MQFix
 * arguments:
 *      1) a pointer to structure MQFix
 *      2) lidx: index of the line
 * return: address to the line
 */
inline uint64_t*
mqfix_line(MQFix* mqf, const uint64_t lidx) {
    assert(MQFIX_FILE == mqf->type);
    return mqf->val + lidx;
}

/* function: mqfix_lvar
 * usage: given index to a line and the index of the target varible, return the
 *      value for the variable
 * arguments:
 *      1) a pointer to structure MQFix
 *      2) lidx: index to the line
 *      3) idx: index of the variable. For x_{n-f+1}, pass 0, and so on
 * return: the value for that variable
 */
inline bool
mqfix_lvar(MQFix* mqf, const uint64_t lidx, const uint64_t idx) {
    uint64_t line = 0x0UL;
    if(MQFIX_FILE == mqf->type) {
        line = *mqfix_line(mqf, lidx);
    } else if(MQFIX_SEQUENTIAL == mqf->type) {
        line = mqf->start + lidx;
    }

    return (line >> idx) & 0x1UL;
}

/* function: mqfix_ltob
 * usage: given index to a line, convert it to an array of bool
 * arguments:
 *      1) a pointer to structure MQFix
 *      2) lidx: index to the line
 *      3) dst: container for the bool array, must be at least as large as the
 *              number of fnum
 * return: void
 */
inline void
mqfix_ltob(MQFix* mqf, const uint64_t lidx, bool* const dst) {
    for(uint64_t i = 0; i < mqf->fnum; ++i) {
        dst[i] = mqfix_lvar(mqf, lidx, i);
    }
}
