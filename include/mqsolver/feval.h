/* feval.h: header file for structure Feval
 */

#ifndef __FEVAL_H__
#define __FEVAL_H__

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>


/* ===========================================================================
 *                              structure Feval definition
 * ===========================================================================
 */
typedef struct {
    uint64_t fnum;  // number of variables to fix
    uint64_t knum;  // number of variables to keep
    uint64_t deg; // degree of the sub-system, same as Macaulay matrix
    uint64_t* d2_offset; // point to the start of 2nd order partial deriv block
    uint64_t* d3_offset; // point to the start of 3rd order partial deriv block
    uint64_t* d4_offset; // point to the start of 4th order partial deriv block
    uint64_t d3_width; // the number of uint64_t used for 3rd order partial deriv
    uint64_t ev[];
    /* ev: an array of uint64_t, which are the evaluation of functions
     *      or partial derivatives of the functions stored in the following
     *      manner:
     *
     *      1. stored column-wise, i.e. an uint64_t holds the term for all the
     *              functions.
     *
     *      2. the first k+1 uint64_t holds the evaluation of f
     *
     *      3. the next fvar_num * (k+1) uint64_t holds the evaluations of 1st
     *              order partial derivatives in grlex order, i.e. x0, x1, ...xf.
     *
     *      4. the next \choose{fvar_num}{2} * (k+1) uint64_t hold the evaluations
     *              of 2nd order partial derivatives in grlex order, i.e. x0x1,
     *              x0x2, x1x2, x0x3, x1x3, x2x3, ..., x0xf, x1xf, ..., x_{f-1}xf
     *
     *      5. the next \choose{fvar_num}{3} * (k+1) uint64_t hold the evaluations
     *              of 3rd partial derivatives in grlex order, i.e. x0x1x2,
     *              x0x1x3, x0x2x3 x1x2x3, ...
     *                      
     *              If deg is 3 then those 3rd partial derivatives are constants
     *              therefore only \choose{fvar_num}{3} uint64_t are used.
     *
     *      6. If deg is 4, then the next \choose{fvar_num}{4} uint64_t hold the
     *              4th order partial derivatives in grlex order
     */
} Feval;

/* ===========================================================================
 *                              function prototypes
 * ===========================================================================
 */

/* function: fev_memsize
 * usage: compute the size of memory needed for the structure
 * arguments:
 *      1) fvar_num: number of variables to fix
 *      2) kvar_num: number of variables to keep
 *      3) deg: degree of sub-system, same as the Macaulay matrix, can be 3, 4
 * return: size of memory needed in bytes
 */
uint64_t
fev_memsize(const uint64_t fvar_num, const uint64_t kvar_num,
            const uint64_t deg);

/* function: fev_create
 * usage: create a Feval container
 * arguments:
 *      1) fvar_num: number of variables to fix
 *      2) kvar_num: number of variables to keep
 *      3) deg: degree of sub-system, same as the Macaulay matrix, can be 3, 4
 * return: a pointer to structure Feval
 */
Feval*
fev_create(const uint64_t fvar_num, const uint64_t kvar_num,
           const uint64_t deg);

/* function: fev_free
 * usage: release a structure Feval
 * arguments:
 *      1) fev: pointer to structure Feval
 * return: void
 */
void
fev_free(Feval* fev);

/* function: fev_getf
 * usage: return the starting address of evaluation of the functions
 * arguments:
 *      1) fev: pointer to structure Feval
 *      2) indices: sorted indices of the vars at the bottom of partial dervative
 *      3) pvar_num: number of vars at the bottom of partial dervative
 * return: the address as uint64_t*
 */
uint64_t*
fev_getf(Feval* fev, const uint64_t* const indices, const uint64_t pvar_num);

/* function: fev_f
 * usage: return the starting address of evaluation of the functions
 * arguments:
 *      1) fev: pointer to structure Feval
 * return the address as uint64_t*
 */
uint64_t*
fev_f(Feval* fev);

/* function: fev_pf1
 * usage: return the starting address of evaluation of the 1st order partial
 *      derivatives of the functions against 1 variable
 * arguments:
 *      1) fev: pointer to structure Feval
 *      2) var_idx: the variable index, should be 0 ~ fvar_num-1
 * return the address as uint64_t*
 */
uint64_t*
fev_pf1(Feval* fev, const uint64_t var_idx);

/* function: fev_pf2
 * usage: return the starting address of evaluation of the 2nd order partial
 *      derivatives of the functions against 2 variables
 * arguments:
 *      1) fev: pointer to structure Feval
 *      2) var1_idx: the smaller variable index, should be 0 ~ fvar_num-2
 *      3) var2_idx: the larger variable index, should be 1 ~ fvar_num-1
 * return the address as uint64_t*
 */
uint64_t*
fev_pf2(Feval* fev, const uint64_t var1_idx, const uint64_t var2_idx);

/* function: fev_pf3
 * usage: return the starting address of evaluation of the 3rd order partial
 *      derivatives of the functions against 3 variables
 * arguments:
 *      1) fev: pointer to structure Feval
 *      2) var1_idx: the smaller variable index, should be 0 ~ fvar_num-3
 *      3) var2_idx: the middle variable index, should be 1 ~ fvar_num-2
 *      4) var3_idx: the largest variable index, should be 2 ~ fvar_num-1
 * return the address as uint64_t*
 */
uint64_t*
fev_pf3(Feval* fev, const uint64_t var1_idx, const uint64_t var2_idx,
        const uint64_t var3_idx);

/* function: fev_pf4
 * usage: return the starting address of evaluation of the 4th order partial
 *      derivatives of the functions against 4 variables. Should only be
 *      called when degree of the sub-system is 4
 * arguments:
 *      1) fev: pointer to structure Feval
 *      2) var1_idx: the smallest variable index, should be 0 ~ fvar_num-4
 *      3) var2_idx: the middle variable index, should be 1 ~ fvar_num-3
 *      4) var3_idx: the 2nd largest variable index, should be 2 ~ fvar_num-2
 *      5) var4_idx: the largest variable index, should be 3 ~ fvar_num-1
 * return the address as uint64_t*
 */
uint64_t*
fev_pf4(Feval* fev, const uint64_t var1_idx, const uint64_t var2_idx,
        const uint64_t var3_idx, const uint64_t var4_idx);

#endif  /* __FEVAL_H__ */
