/* fix.h: header file that defines functions for fixing variables
 */

#ifndef __FIX_H__
#define __FIX_H__

#include <stdint.h>
#include <assert.h>

#include "cuda_util.h"
#include "graycode.h"

/* macro: eterm_lsys
 * usage: given the thread id in the warp, the interleaved storage for lsys,
 *      and the term index, extract the term
 * arguments:
 *      1) lsys: interleaved storage for the warp
 *      2) warp_tid: thread id in the warp
 *      3) idx: index specifying the term, for x1, pass 0, enz. For constant
 *              term, pass the number of variables to keep
 * return: void
 */
#define eterm_lsys(lsys, warp_tid, idx) \
    ((lsys)[(idx) * WARP_SIZE + (warp_tid)])

/* macro: add_lsys
 * usage: given the thread id in the warp, the interleaved storage for lsys,
 *      and the term to add, add the term to the right slot
 * arguments:
 *      1) lsys: interleaved storage for the warp
 *      2) warp_tid: thread id in the warp
 *      3) idx: index specifying the term, for x1, pass 0, enz. For constant
 *              term, pass the number of variables to keep
 *      4) term: term to fill
 * return: void
 */
#define add_lsys(lsys, warp_tid, idx, term) do { \
    eterm_lsys(lsys, warp_tid, idx) ^= (term); \
} while(0)

/* function: fix_pf[1-2]
 * usage: fix variables in the partial derivative of a subsytem and turn it
 *      into a linear system. Note that the subsytem must be transposed and
 *      stored as an array of uint32_t where each uint32_t represents the same
 *      column for all rows. The max number of eqs is therefore 32.
 *      The non-linear terms should be kept in the subsystem even though they
 *      are all zeros.
 * arguments:
 *      1) lsys: container for the linear system for the warp, transposed,
 *              interleaved and should be clear to zero
 *      2) sub: the sub-system
 *      3) deg: degree of the sub-system, same as the original Macaulay matrix
 *      4) vnum: number of variables
 *      5) knum: number of variables to keep
 *      6) kfnum: number of variables to fix by the kernel
 *      7) tv: the values of the variables to fix by a thread, stored as a
 *              32-bit integer where the last bit is the value for x_{n-kf}
 *      8) kv: same as tv but for kf vars. The last bit is the value for x_n
 *      9+) idx[1-2]?: sorted indices of the variables w.r.t. which the paritial
 *              derivatives are taken.
 */

/* fix_pf1: fix 1st order paritial derivatives */
__device__ void
fix_pf1(uint32_t* const __restrict__ lsys,
        const uint32_t* const __restrict__ sub, const uint32_t deg,
        const uint32_t vnum, const uint32_t knum, const uint32_t kfnum,
        const uint32_t tv, const uint32_t kv, const uint32_t idx);

/* fix_pf2: fix 2nd order paritial derivatives */
__device__ void
fix_pf2(uint32_t* const __restrict__ lsys,
        const uint32_t* const __restrict__ sub, const uint32_t deg,
        const uint32_t vnum, const uint32_t knum, const uint32_t kfnum,
        const uint32_t tv, const uint32_t kv, const uint32_t idx1,
        const uint32_t idx2);

/* function: fix_pf3
 * usage: fix variables in the partial derivative of a subsytem and turn it
 *      into a linear system. Note that the subsytem must be transposed and
 *      stored as an array of uint32_t where each uint32_t represents the same
 *      column for all rows. The max number of eqs is therefore 32.
 *      The non-linear terms should be kept in the subsystem even though they
 *      are all zeros. Only for Macaulay degree 4 since degree 3 is trivial.
 * arguments:
 *      1) lsys: container for the linear system for the warp, transposed,
 *              interleaved and should be clear to zero
 *      2) sub: the sub-system
 *      3) vnum: number of variables
 *      4) knum: number of variables to keep
 *      5) kfnum: number of variables to fix by the kernel
 *      6) tv: the values of the variables to fix by a thread, stored as a
 *              32-bit integer where the last bit is the value for x_{n-kf}
 *      7) kv: same as tv but for kf vars. The last bit is the value for x_n
 *      8+) idx[1-3]: sorted indices of the variables w.r.t. which the paritial
 *              derivatives are taken.
 */
__device__ void
fix_pf3(uint32_t* const __restrict__ lsys,
        const uint32_t* const __restrict__ sub,
        const uint32_t vnum, const uint32_t knum, const uint32_t kfnum,
        const uint32_t tv, const uint32_t kv, const uint32_t idx1,
        const uint32_t idx2, const uint32_t idx3);

/* function: fix_f
 * usage: fix variables in a subsystem for graycode enumeration. Note
 *      that the subsystem must be transposed and stored as an array of
 *      uint32_t where each uint32_t represents the same column for all rows.
 *      The max number of eqs is therefore 32. The non-linear terms should be
 *      kept in the sub-system. Out of all variables to fix, last kfnum
 *      variables are set according to v, while the rest are all set to 0.
 * arguments:
 *      1) lsys: container for the linear system for the warp, transposed and
 *              interleaved and should be clear to zero
 *      2) sub: the sub-system
 *      3) deg: degree of the sub-system, same as the original Macaulay matrix
 *      4) vnum: number of variables
 *      5) knum: number of variables to keep
 *      6) kfnum: number of variables to fix by the kernel
 *      7) v: the values of the variables to fix by the kernel, stored as a 32-bit
 *              integer where the last bit is the value for x_n
 * return: void
 */
__device__ void
fix_f(uint32_t* const __restrict__ lsys, const uint32_t* const __restrict__ sub,
      const uint32_t deg, const uint32_t vnum, const uint32_t knum,
      const uint32_t kfnum, const uint32_t v);

/* function: fix_sub
 * usage: fix variables in a subsystem to reduce the number of variables.
 *      The subsystem must be transposed and stored as an array of
 *      uint32_t where each uint32_t represents the same column for all rows.
 *      The max number of eqs is therefore 32. The non-linear terms should be
 *      kept in the sub-system.
 * arguments:
 *      1) dst: container for the resultant sub-system, should be large enough
 *              to hold the number of remaining monomials
 *      2) sub: the sub-system
 *      3) deg: degree of the sub-system, can be 3 or 4
 *      4) vnum: number of variables in sub
 *      5) fnum: number of variables to fix, starting from the last one
 *      6) v: the values of the variables to fix by the kernel, stored as a 32-bit
 *              integer where the last bit is the value for x_n
 * return: void
 */
__host__ void
fix_sub(uint32_t* const __restrict__ dst, const uint32_t* const __restrict__ sub,
        const uint32_t deg, const uint32_t vnum, const uint32_t fnum, const uint32_t v);

#endif  /* __FIX_H__ */
