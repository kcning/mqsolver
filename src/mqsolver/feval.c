/* feval.c: implementation of structure Feval
 */

#include "feval.h"

/* subroutines */
static inline uint64_t
binom2(const uint64_t n) {
    return n * (n-1) / 2;
}

static inline uint64_t
binom3(const uint64_t n) {
    return n * (n-1) * (n-2) / 6;
}

static inline uint64_t
binom4(const uint64_t n) {
    return n * (n-1) * (n-2) * (n-3) / 24;
}

/* function: fev_slot_num
 * usage: compute the number of uint64_t needed for the structure
 * arguments:
 *      1) fvar_num: number of variables to fix
 *      2) kvar_num: number of variables to keep
 *      3) deg: degree of sub-system, same as the Macaulay matrix, can be 3, 4
 * return: number of uint64_t needed
 */
static inline uint64_t
fev_slot_num(const uint64_t fvar_num, const uint64_t kvar_num,
             const uint64_t deg) {
    assert(3 == deg || 4 == deg);

    // compute the number of uint64_t needed
    uint64_t int_num = (kvar_num+1); // f
    int_num += fvar_num * (kvar_num+1); // 1st order partial dervs
    int_num += binom2(fvar_num) * (kvar_num+1); // 2nd order

    if(3 == deg) {
        int_num += binom3(fvar_num); // 3rd order, constant
    } else {
        int_num += binom3(fvar_num) * (kvar_num+1); // 3rd order
        int_num += binom4(fvar_num); // 4th order, constant
    }

    return int_num;
}

/* function: fev_memsize
 * usage: compute the size of memory needed for the structure
 * arguments:
 *      1) fvar_num: number of variables to fix
 *      2) kvar_num: number of variables to keep
 *      3) deg: degree of sub-system, same as the Macaulay matrix, can be 3, 4
 * return: size of memory needed in bytes
 */
inline uint64_t
fev_memsize(const uint64_t fvar_num, const uint64_t kvar_num,
            const uint64_t deg) {
    return fev_slot_num(fvar_num, kvar_num, deg) * sizeof(uint64_t) + sizeof(Feval);
}


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
           const uint64_t deg) {
    assert(3 == deg || 4 == deg);
    const uint64_t int_num = fev_slot_num(fvar_num, kvar_num, deg);

    Feval* fev = (Feval*) malloc(sizeof(Feval) + sizeof(uint64_t) * int_num);
    if(NULL == fev) {
        return NULL;
    }
    

    fev->fnum = fvar_num;
    fev->knum = kvar_num;
    fev->deg = deg;

    fev->d2_offset = fev->ev + (1 + fev->fnum) * (fev->knum+1);
    fev->d3_offset = fev->d2_offset + binom2(fev->fnum) * (fev->knum+1);

    if(4 == deg) {
        fev->d3_width = fev->knum + 1;
        fev->d4_offset = fev->d3_offset + binom3(fev->fnum) * fev->d3_width;
    } else {
        fev->d3_width = 1;
        fev->d4_offset = NULL;
    }

    return fev;
}

/* function: fev_free
 * usage: release a structure Feval
 * arguments:
 *      1) fev: pointer to structure Feval
 * return: void
 */
inline void
fev_free(Feval* fev) {
    free(fev);
}

/* function: fev_getf
 * usage: return the starting address of evaluation of the functions
 * arguments:
 *      1) fev: pointer to structure Feval
 *      2) indices: sorted indices of the vars at the bottom of partial dervative
 *      3) pvar_num: number of vars at the bottom of partial dervative
 * return: the address as uint64_t*
 */
inline uint64_t*
fev_getf(Feval* fev, const uint64_t* const indices, const uint64_t pvar_num) {
    assert(pvar_num <= fev->deg);

    switch(pvar_num) {
        case 0:
            return fev_f(fev);
        case 1:
            return fev_pf1(fev, indices[0]);
        case 2:
            return fev_pf2(fev, indices[0], indices[1]);
        case 3:
            return fev_pf3(fev, indices[0], indices[1], indices[2]);
        case 4:
            return fev_pf4(fev, indices[0], indices[1], indices[2], indices[3]);
        default:
            return NULL;
    }
}

/* function: fev_f
 * usage: return the starting address of evaluation of the functions
 * arguments:
 *      1) fev: pointer to structure Feval
 * return: the address as uint64_t*
 */
inline uint64_t*
fev_f(Feval* fev) {
    return fev->ev;
}

/* function: fev_pf1
 * usage: return the starting address of evaluation of the 1st order partial
 *      derivatives of the functions against 1 variable
 * arguments:
 *      1) fev: pointer to structure Feval
 *      2) var_idx: the variable index, should be 0 ~ fvar_num-1
 * return: the address as uint64_t*
 */
inline uint64_t*
fev_pf1(Feval* fev, const uint64_t var_idx) {
    return fev->ev + ((1 + var_idx) * (fev->knum + 1));
}

/* function: fev_pf2
 * usage: return the starting address of evaluation of the 2nd order partial
 *      derivatives of the functions against 2 variables
 * arguments:
 *      1) fev: pointer to structure Feval
 *      2) var1_idx: the smaller variable index, should be 0 ~ fvar_num-2
 *      3) var2_idx: the larger variable index, should be 1 ~ fvar_num-1
 * return: the address as uint64_t*
 */
inline uint64_t*
fev_pf2(Feval* fev, const uint64_t var1_idx, const uint64_t var2_idx) {
    uint64_t offset = (var2_idx * (var2_idx-1) / 2 + var1_idx) * (fev->knum+1);
    return fev->d2_offset + offset;
}

/* function: fev_pf3
 * usage: return the starting address of evaluation of the 3rd order partial
 *      derivatives of the functions against 3 variables
 * arguments:
 *      1) fev: pointer to structure Feval
 *      2) var1_idx: the smaller variable index, should be 0 ~ fvar_num-3
 *      3) var2_idx: the middle variable index, should be 1 ~ fvar_num-2
 *      4) var3_idx: the largest variable index, should be 2 ~ fvar_num-1
 * return: the address as uint64_t*
 */
inline uint64_t*
fev_pf3(Feval* fev, const uint64_t var1_idx, const uint64_t var2_idx,
        const uint64_t var3_idx) {
    const uint64_t offset = var3_idx * (var3_idx-1) * (var3_idx-2) / 6
                            + var2_idx * (var2_idx-1) / 2 
                            + var1_idx;
    return fev->d3_offset + offset * fev->d3_width;
}

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
 * return: the address as uint64_t*
 */
inline uint64_t*
fev_pf4(Feval* fev, const uint64_t var1_idx, const uint64_t var2_idx,
        const uint64_t var3_idx, const uint64_t var4_idx) {
    assert(4 == fev->deg);

    const uint64_t offset = var4_idx * (var4_idx-1) * (var4_idx-2) * (var4_idx-3) / 24
                            + var3_idx * (var3_idx-1) * (var3_idx-2) / 6
                            + var2_idx * (var2_idx-1) / 2
                            + var1_idx;
    return fev->d4_offset + offset;
}
