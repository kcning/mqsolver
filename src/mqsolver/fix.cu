/* fix.cu: implementation of fix.h
 */

extern "C" {

#include "fix.h"

/* vprod[1-4]?: subroutines for fix_p?f[1-3]?, given the values of the variables,
 * compute their product.
 *
 * The indices of the variables to fix DOES NOT need to be sorted, e.g. for x4
 * and x6, pass either [3, 5] or [5, 3] is fine. The max size of indices is 4
 * because we only work with Macaulay matrix of degree 3 or 4.
 *
 * The format of v is the same as fix_pf1, note that we can this way kernel
 * and thread can fix at most 32 variables respectively.
 */

__device__ __host__ static __forceinline__ uint32_t
vprod(const uint32_t v, const uint32_t vnum, const uint32_t pvnum,
      const uint32_t idx1, const uint32_t idx2, const uint32_t idx3,
      const uint32_t idx4) {
    uint32_t c1 = (v >> (vnum-1-idx1)) | (pvnum < 1);
    uint32_t c2 = (v >> (vnum-1-idx2)) | (pvnum < 2);
    uint32_t c3 = (v >> (vnum-1-idx3)) | (pvnum < 3);
    uint32_t c4 = (v >> (vnum-1-idx4)) | (pvnum < 4);
    return c1 & c2 & c3 & c4 & 0x1U;
}

__device__ __host__ static __forceinline__ uint32_t
vprod1(const uint32_t v, const uint32_t vnum, const uint32_t idx) {
    return (v >> (vnum-1-idx) & 0x1U);
}

__device__ __host__ static __forceinline__ uint32_t
vprod2(const uint32_t v, const uint32_t vnum,
       const uint32_t idx1, const uint32_t idx2) {
    return (v >> (vnum-1-idx1)) & (v >> (vnum-1-idx2)) & 0x1U;
}

__device__ __host__ static __forceinline__ uint32_t
vprod3(const uint32_t v, const uint32_t vnum, const uint32_t idx1,
       const uint32_t idx2, const uint32_t idx3) {
    uint32_t c1 = v >> (vnum-1-idx1);
    uint32_t c2 = v >> (vnum-1-idx2);
    uint32_t c3 = v >> (vnum-1-idx3);
    return c1 & c2 & c3 & 0x1U;
}

__device__ __host__ static __forceinline__ uint32_t
vprod4(const uint32_t v, const uint32_t vnum, const uint32_t idx1,
       const uint32_t idx2, const uint32_t idx3, const uint32_t idx4) {
    uint32_t c1 = v >> (vnum-1-idx1);
    uint32_t c2 = v >> (vnum-1-idx2);
    uint32_t c3 = v >> (vnum-1-idx3);
    uint32_t c4 = v >> (vnum-1-idx4);
    return c1 & c2 & c3 & c4 & 0x1U;
}

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
        const uint32_t tv, const uint32_t kv, const uint32_t idx) {

    const uint32_t front_num = vnum - kfnum;
    const uint32_t warp_tid = warp_tid();

    if(4 == deg) { // degree 4 monomials
        // fix 1st var
        // reduce to constant
        // last var fix by thread
        for(uint32_t i = idx+3; i < front_num; ++i) {

            if(!vprod1(tv, front_num, i)) {
                continue;
            }

            // fix by thread
            for(uint32_t j = idx+2; j < i; ++j) {
                
                if(!vprod1(tv, front_num, j)) {
                    continue;
                }

                for(uint32_t k = idx+1; k < j; ++k) {
                    if(vprod1(tv, front_num, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx, k, j, i)]);
                    }
                }
            }
        }

        // last var fix by kernel
        for(uint32_t i = front_num; i < vnum; ++i) {

            if(!vprod1(kv, vnum, i)) {
                continue;
            }

            // 2nd last var fix by thread
            for(uint32_t j = idx+2; j < front_num; ++j) {

                if(!vprod1(tv, front_num, j)) {
                    continue;
                }

                for(uint32_t k = idx+1; k < j; ++k) {
                    if(vprod1(tv, front_num, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx, k, j, i)]);
                    }
                }
            }

            // 2nd last var fix by kernel
            for(uint32_t j = front_num; j < i; ++j) {

                if(!vprod1(kv, vnum, j)) {
                    continue;
                }
                    
                // 2nd var fix by thread
                for(uint32_t k = idx+1; k < front_num; ++k) {
                    if(vprod1(tv, front_num, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx, k, j, i)]);
                    }
                }

                // 2nd var fix by kernel
                for(uint32_t k = front_num; k < j; ++k) {
                    if(vprod1(kv, vnum, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx, k, j, i)]);
                    }
                }
            }
        }

        // fix 2nd var
        // last var fix by thread
        for(uint32_t i = idx+2; i < front_num; ++i) {

            if(!vprod1(tv, front_num, i)) {
                continue;
            }

            // fix by thread
            for(uint32_t j = idx+1; j < i; ++j) {
            
                if(!vprod1(tv, front_num, j)) {
                    continue;
                }

                // reduce to linear
                for(uint32_t k = 0; k < knum; ++k) {
                    add_lsys(lsys, warp_tid, k, sub[midx4(vnum, k, idx, j, i)]);
                }

                // reduce to constant
                for(uint32_t k = knum; k < idx; ++k) {
                    if(vprod1(tv, front_num, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, k, idx, j, i)]);
                    }
                }
            }
        }

        // last var fix by kernel
        for(uint32_t i = front_num; i < vnum; ++i) {

            if(!vprod1(kv, vnum, i)) {
                continue;
            }

            // 2nd last var fix by thread
            for(uint32_t j = idx+1; j < front_num; ++j) {
            
                if(!vprod1(tv, front_num, j)) {
                    continue;
                }

                // reduce to linear
                for(uint32_t k = 0; k < knum; ++k) {
                    add_lsys(lsys, warp_tid, k, sub[midx4(vnum, k, idx, j, i)]);
                }

                // reduce to constant
                for(uint32_t k = knum; k < idx; ++k) {
                    if(vprod1(tv, front_num, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, k, idx, j, i)]);
                    }
                }
            }

            // 2nd last fix by kernel
            for(uint32_t j = front_num; j < i; ++j) {

                if(!vprod1(kv, vnum, j)) {
                    continue;
                }

                // reduce to linear
                for(uint32_t k = 0; k < knum; ++k) {
                    add_lsys(lsys, warp_tid, k, sub[midx4(vnum, k, idx, j, i)]);
                }

                // reduce to constant
                for(uint32_t k = knum; k < idx; ++k) {
                    if(vprod1(tv, front_num, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, k, idx, j, i)]);
                    }
                }
            }
        }

        // fix 3rd var
        // last var fix by thread
        for(uint32_t i = idx+1; i < front_num; ++i) {

            if(!vprod1(tv, front_num, i)) {
                continue;
            }

            for(uint32_t j = knum; j < idx; ++j) {
                
                if(!vprod1(tv, front_num, j)) {
                    continue;
                }

                // reduce to linear
                for(uint32_t k = 0; k < knum; ++k) {
                    add_lsys(lsys, warp_tid, k, sub[midx4(vnum, k, j, idx, i)]);
                }

                // reduce to constant
                for(uint32_t k = knum; k < j; ++k) {
                    if(vprod1(tv, front_num, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, k, j, idx, i)]);
                    }
                }
            }
        }

        // last var fix by kernel
        for(uint32_t i = front_num; i < vnum; ++i) {

            if(!vprod1(kv, vnum, i)) {
                continue;
            }

            // fix by thread
            for(uint32_t j = knum; j < idx; ++j) {
                
                if(!vprod1(tv, front_num, j)) {
                    continue;
                }

                // reduce to linear
                for(uint32_t k = 0; k < knum; ++k) {
                    add_lsys(lsys, warp_tid, k, sub[midx4(vnum, k, j, idx, i)]);
                }

                // reduce to constant
                for(uint32_t k = knum; k < j; ++k) {
                    if(vprod1(tv, front_num, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, k, j, idx, i)]);
                    }
                }
            }
        }

        // fix 4th var, all vars fix by thread
        for(uint32_t i = knum; i < idx; ++i) {

            if(!vprod1(tv, front_num, i)) {
                continue;
            }

            for(uint32_t j = knum; j < i; ++j) {

                if(!vprod1(tv, front_num, j)) {
                    continue;
                }

                // reduce to linear
                for(uint32_t k = 0; k < knum; ++k) {
                    add_lsys(lsys, warp_tid, k, sub[midx4(vnum, k, j, i, idx)]);
                }

                // reduce to constant
                for(uint32_t k = knum; k < j; ++k) {
                    if(vprod1(tv, front_num, k)) {
                        add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, k, j, i, idx)]);
                    }
                }
            }
        }
    } // end of degree 4 monomials

    // degree 3 monomials
    // fix 1st var
    // reduce to constant
    // last var fixed by thread
    for(uint32_t i = idx+2; i < front_num; ++i) {

        if(!vprod1(tv, front_num, i)) {
            continue;
        }

        for(uint32_t j = idx+1; j < i; ++j) {
            if(vprod1(tv, front_num, j)) {
                add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, idx, j, i)]);
            }
        }
    }

    // last var fixed by kernel
    for(uint32_t i = front_num; i < vnum; ++i) {

        if(!vprod1(kv, vnum, i)) {
            continue;
        }

        // fix by thread
        for(uint32_t j = idx+1; j < front_num; ++j) {
            if(vprod1(tv, front_num, j)) {
                add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, idx, j, i)]);
            }
        }

        // fix by kernel
        for(uint32_t j = front_num; j < i; ++j) {
            if(vprod1(kv, vnum, j)) {
                add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, idx, j, i)]);
            }
        }
    }

    // fix 2nd var
    // last var fix by thread
    for(uint32_t i = idx+1; i < front_num; ++i) {

        if(!vprod1(tv, front_num, i)) {
            continue;
        }

        // reduce to linear
        for(uint32_t j = 0; j < knum; ++j) {
            add_lsys(lsys, warp_tid, j, sub[midx3(deg, vnum, j, idx, i)]);
        }

        // reduce to constant
        for(uint32_t j = knum; j < idx; ++j) {
            if(vprod1(tv, front_num, j)) {
                add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, j, idx, i)]);
            }
        }
    }

    // last var fix by kernel
    for(uint32_t i = front_num; i < vnum; ++i) {

        if(!vprod1(kv, vnum, i)) {
            continue;
        }

        // reduce to linear
        for(uint32_t j = 0; j < knum; ++j) {
            add_lsys(lsys, warp_tid, j, sub[midx3(deg, vnum, j, idx, i)]);
        }

        // reduce to constant
        for(uint32_t j = knum; j < idx; ++j) {
            if(vprod1(tv, front_num, j)) {
                add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, j, idx, i)]);
            }
        }
    }

    // fix 3rd var
    // fix by thread
    for(uint32_t i = knum; i < idx; ++i) {

        if(!vprod1(tv, front_num, i)) {
            continue;
        }

        // reduce to linear 
        for(uint32_t j = 0; j < knum; ++j) {
            add_lsys(lsys, warp_tid, j, sub[midx3(deg, vnum, j, i, idx)]);
        }

        // reduce to constant
        for(uint32_t j = knum; j < i; ++j) {
            if(vprod1(tv, front_num, j)) {
                add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, j, i, idx)]);
            }
        }
    }

    // degree 2 monomials
    // free var is 1st
    // reduce to linear
    for(uint32_t i = 0; i < knum; ++i) {
        add_lsys(lsys, warp_tid, i, sub[midx2(deg, vnum, i, idx)]);
    }

    // reduce to constant
    for(uint32_t i = knum; i < idx; ++i) {
        if(vprod1(tv, front_num, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx2(deg, vnum, i, idx)]);
        }
    }

    // free var is 2nd
    // reduce to constant
    // fix by thread
    for(uint32_t i = idx+1; i < front_num; ++i) {
        if(vprod1(tv, front_num, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx2(deg, vnum, idx, i)]);
        }
    }

    // fix by kernel
    for(uint32_t i = front_num; i < vnum; ++i) {
        if(vprod1(kv, vnum, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx2(deg, vnum, idx, i)]);
        }
    }
    
    // degree 1 monomials
    add_lsys(lsys, warp_tid, knum, sub[midx1(deg, vnum, idx)]);
}

/* fix_pf2: fix 2nd order paritial derivatives */
__device__ void
fix_pf2(uint32_t* const __restrict__ lsys,
        const uint32_t* const __restrict__ sub, const uint32_t deg,
        const uint32_t vnum, const uint32_t knum, const uint32_t kfnum,
        const uint32_t tv, const uint32_t kv, const uint32_t idx1,
        const uint32_t idx2) {

    const uint32_t front_num = vnum - kfnum;
    const uint32_t warp_tid = warp_tid();

    if(4 == deg) { // degree 4 monomials
        // free vars are 1st and 2nd
        // vars fix by thread
        for(uint32_t i = knum; i < idx1; ++i) {

            if(!vprod1(tv, front_num, i)) {
                continue;
            }

            // reduce to linear term
            for(uint32_t j = 0; j < knum; ++j) {
                add_lsys(lsys, warp_tid, j, sub[midx4(vnum, j, i, idx1, idx2)]);
            }

            // reduce to constant
            for(uint32_t j = knum; j < i; ++j) {
                if(vprod1(tv, front_num, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, j, i, idx1, idx2)]);
                }
            }
        }

        // free vars are 1st and 3rd
        // fix by thread
        for(uint32_t i = idx1+1; i < idx2; ++i) {

            if(!vprod1(tv, front_num, i)) {
                continue;
            }

            // reduce to linear term
            for(uint32_t j = 0; j < knum; ++j) {
                add_lsys(lsys, warp_tid, j, sub[midx4(vnum, j, idx1, i, idx2)]);
            }

            // reduce to constant
            for(uint32_t j = knum; j < idx1; ++j) {
                if(vprod1(tv, front_num, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, j, idx1, i, idx2)]);
                }
            }
        }

        // free vars are 1st and 4th
        // fix by thread
        for(uint32_t i = idx2+1; i < front_num; ++i) {

            if(!vprod1(tv, front_num, i)) {
                continue;
            }
            
            // reduce to linear term
            for(uint32_t j = 0; j < knum; ++j) {
                add_lsys(lsys, warp_tid, j, sub[midx4(vnum, j, idx1, idx2, i)]);
            }

            // reduce to constant
            for(uint32_t j = knum; j < idx1; ++j) {
                if(vprod1(tv, front_num, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, j, idx1, idx2, i)]);
                }
            }
        }

        // fix by kernel
        for(uint32_t i = front_num; i < vnum; ++i) {
            
            if(!vprod1(kv, vnum, i)) {
                continue;
            }

            // reduce to linear term
            for(uint32_t j = 0; j < knum; ++j) {
                add_lsys(lsys, warp_tid, j, sub[midx4(vnum, j, idx1, idx2, i)]);
            }

            // reduce to constant
            for(uint32_t j = knum; j < idx1; ++j) {
                if(vprod1(tv, front_num, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, j, idx1, idx2, i)]);
                }
            }
        }

        // free vars are 2nd and 3rd
        // reduce to constant
        // fix by thread
        for(uint32_t i = idx1+2; i < idx2; ++i) {

            if(!vprod1(tv, front_num, i)) {
                continue;
            }

            for(uint32_t j = idx1+1; j < i; ++j) {
                if(vprod1(tv, front_num, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, j, i, idx2)]);
                }
            }
        }

        // free vars are 2nd and 4th
        // reduce to constant
        // fix by thread
        for(uint32_t i = idx2+1; i < front_num; ++i) {
            if(!vprod1(tv, front_num, i)) {
                continue;
            }

            for(uint32_t j = idx1+1; j < idx2; ++j) {
                if(vprod1(tv, front_num, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, j, idx2, i)]);
                }
            }
        }

        // fix by kernel
        for(uint32_t i = front_num; i < vnum; ++i) {
            if(!vprod1(kv, vnum, i)) {
                continue;
            }

            for(uint32_t j = idx1+1; j < idx2; ++j) {
                if(vprod1(tv, front_num, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, j, idx2, i)]);
                }
            }
        }

        // free vars are 3rd and 4th
        // reduce to constant
        // fix by thread
        for(uint32_t i = idx2+2; i < front_num; ++i) {

            if(!vprod1(tv, front_num, i)) {
                continue;
            }

            for(uint32_t j = idx2+1; j < i; ++j) {
                if(vprod1(tv, front_num, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, idx2, j, i)]);
                }
            }
        }

        // fix by kernel
        for(uint32_t i = front_num; i < vnum; ++i) {

            if(!vprod1(kv, vnum, i)) {
                continue;
            }

            // fix by thread
            for(uint32_t j = idx2+1; j < front_num; ++j) {
                if(vprod1(tv, front_num, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, idx2, j, i)]);
                }
            }

            // fix by kernel
            for(uint32_t j = front_num; j < i; ++j) {
                if(vprod1(kv, vnum, j)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, idx2, j, i)]);
                }
            }
        }
    } // end of degree 4 monomials

    // degree 3 monomials
    
    // free var is 1st one
    // reduce to linear term
    // fix by thread
    for(uint32_t i = 0; i < knum; ++i) {
        add_lsys(lsys, warp_tid, i, sub[midx3(deg, vnum, i, idx1, idx2)]);
    }

    // reduce to constant
    for(uint32_t i = knum; i < idx1; ++i) {
        if(vprod1(tv, front_num, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, i, idx1, idx2)]);
        }
    }

    // free var is 2nd one
    // reduce to constant
    // fix by thread
    for(uint32_t i = idx1+1; i < idx2; ++i) {
        if(vprod1(tv, front_num, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, idx1, i, idx2)]);
        }
    }

    // free var is 3rd one
    // reduce to constant
    // fix by thread
    for(uint32_t i = idx2+1; i < front_num; ++i) {
        if(vprod1(tv, front_num, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, idx1, idx2, i)]);
        }
    }

    // fix by kernel
    for(uint32_t i = front_num; i < vnum; ++i) {
        if(vprod1(kv, vnum, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, idx1, idx2, i)]);
        }
    }

    // degree 2 monomials
    // reduce to constant
    add_lsys(lsys, warp_tid, knum, sub[midx2(deg, vnum, idx1, idx2)]);
}

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
 *      8+) vars: sorted indices of the variables w.r.t. which the paritial
 *              derivatives are taken.
 */
__device__ void
fix_pf3(uint32_t* const __restrict__ lsys,
        const uint32_t* const __restrict__ sub,
        const uint32_t vnum, const uint32_t knum, const uint32_t kfnum,
        const uint32_t tv, const uint32_t kv, const uint32_t idx1,
        const uint32_t idx2, const uint32_t idx3) {

    const uint32_t front_num = vnum - kfnum;
    const uint32_t warp_tid = warp_tid();
    // degree 4 monomials

    // free var is 1st one
    // reduce to linear term
    for(uint32_t i = 0; i < knum; ++i) {
        add_lsys(lsys, warp_tid, i, sub[midx4(vnum, i, idx1, idx2, idx3)]);
    }

    // reduce to constant
    for(uint32_t i = knum; i < idx1; ++i) {
        if(vprod1(tv, front_num, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, i, idx1, idx2, idx3)]);
        }
    }

    // free var is 2nd one
    // reduce to constant
    for(uint32_t i = idx1+1; i < idx2; ++i) {
        if(vprod1(tv, front_num, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, i, idx2, idx3)]);
        }
    }

    // free var is 3rd one
    // reduce to constant
    for(uint32_t i = idx2+1; i < idx3; ++i) {
        if(vprod1(tv, front_num, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, idx2, i, idx3)]);
        }
    }

    // free var is 4th one
    // reduce to constant
    // vars fixed by thread
    for(uint32_t i = idx3+1; i < front_num; ++i) {
        if(vprod1(tv, front_num, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, idx2, idx3, i)]);
        }
    }

    // vars fixed by kernel
    // reduce to constant
    for(uint32_t i = front_num; i < vnum; ++i) {
        if(vprod1(kv, vnum, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, idx1, idx2, idx3, i)]);
        }
    }

    // degree 3 monomials
    add_lsys(lsys, warp_tid, knum, sub[midx3(4, vnum, idx1, idx2, idx3)]);
}

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
      const uint32_t kfnum, const uint32_t v) {

    const uint32_t front_num = vnum - kfnum; // num of vars before last kfnum vars
    const uint32_t flterm_idx = midx1(deg, vnum, 0); // index to x1
    const uint32_t warp_tid = warp_tid();

    // copy first k linear terms and constant term
    for(uint32_t i = 0; i < knum; ++i) {
        add_lsys(lsys, warp_tid, i, sub[flterm_idx+i]);
    }
    add_lsys(lsys, warp_tid, knum, sub[flterm_idx + vnum]);

    if(4 == deg) {
        for(int i = front_num+2; i < vnum; ++i) {
            for(int j = front_num+1; j < i; ++j) {
                for(int k = front_num; k < j; ++k) {
                    for(int l = 0; l < knum; ++l) {
                        // become linear
                        if(vprod3(v, vnum, i, j, k)) {
                            add_lsys(lsys, warp_tid, l, sub[midx4(vnum, l, k, j, i)]);
                        }
                    }

                    for(int l = front_num; l < k; ++l) {
                        // become constant
                        if(vprod4(v, vnum, i, j, k, l)) {
                            add_lsys(lsys, warp_tid, knum, sub[midx4(vnum, l, k, j, i)]);
                        }
                    }
                }
            }
        }
    }

    // degree 3
    for(int i = front_num+1; i < vnum; ++i) {
        for(int j = front_num; j < i; ++j) {
            for(int k = 0; k < knum; ++k) {
                // become linear
                if(vprod2(v, vnum, i, j)) {
                    add_lsys(lsys, warp_tid, k, sub[midx3(deg, vnum, k, j, i)]);
                }
            }

            for(int k = front_num; k < j; ++k) {
                if(vprod3(v, vnum, i, j, k)) {
                    add_lsys(lsys, warp_tid, knum, sub[midx3(deg, vnum, k, j, i)]);
                }
            }
        }
    }

    // degree 2
    for(int i = front_num; i < vnum; ++i) {
        for(int j = 0; j < knum; ++j) {
            // become linear
            if(vprod1(v, vnum, i)) {
                add_lsys(lsys, warp_tid, j, sub[midx2(deg, vnum, j, i)]);
            }
        }

        for(int j = front_num; j < i; ++j) {
            // become constant
            if(vprod2(v, vnum, i, j)) {
                add_lsys(lsys, warp_tid, knum, sub[midx2(deg, vnum, j, i)]);
            }
        }
    }
    
    // degree 1
    for(int i = front_num; i < vnum; ++i) {
        // become constant
        if(vprod1(v, vnum, i)) {
            add_lsys(lsys, warp_tid, knum, sub[midx1(deg, vnum, i)]);
        }
    }
}

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
        const uint32_t deg, const uint32_t vnum, const uint32_t fnum, const uint32_t v) {

    if(!fnum) { // fix no variables
        memcpy(dst, sub, sizeof(uint32_t) * (midx0(deg, vnum) + 1));
        return;
    }

    uint32_t rvnum = vnum - fnum;
    uint32_t ltnum = midx0(deg, rvnum) + 1;
    memset(dst, 0x0, sizeof(uint32_t) * ltnum);

    if(4 == deg) {
        ltnum = binom4(rvnum);
        memcpy(dst, sub, sizeof(uint32_t) * ltnum);

        // the rest
        for(uint32_t i = rvnum; i < vnum; ++i) {
            if(!vprod1(v, vnum, i)) {
                continue;
            }

            // reduce to degree 3 monomials
            for(uint32_t j = 2; j < rvnum; ++j) {
                for(uint32_t k = 1; k < j; ++k) {
                    for(uint32_t l = 0; l < k; ++l) {
                        dst[midx3(4, rvnum, l, k, j)] ^= sub[midx4(vnum, l, k, j, i)];
                    }
                }
            }

            for(uint32_t j = rvnum; j < i; ++j) {
                if(!vprod1(v, vnum, j)) {
                    continue;
                }

                // reduce to degree 2 monomials
                for(uint32_t k = 1; k < rvnum; ++k) {
                    for(uint32_t l = 0; l < k; ++l) {
                        dst[midx2(4, rvnum, l, k)] ^= sub[midx4(vnum, l, k, j, i)];
                    }
                }

                for(uint32_t k = rvnum; k < j; ++k) {
                    if(!vprod1(v, vnum, k)) {
                        continue;
                    }

                    // reduce to linear term
                    for(uint32_t l = 0; l < rvnum; ++l) {
                        dst[midx1(4, rvnum, l)] ^= sub[midx4(vnum, l, k, j, i)];
                    }

                    // reduce to constant term
                    for(uint32_t l = rvnum; l < k; ++l) {
                        if(vprod1(v, vnum, l)) {
                            dst[midx0(4, rvnum)] ^= sub[midx4(vnum, l, k, j, i)];
                        }
                    }
                }
            }
        }
    }

    // degree 3 monomials
    uint32_t dst_offset = midx3(deg, rvnum, 0, 1, 2);
    uint32_t src_offset = midx3(deg, vnum, 0, 1, 2);
    ltnum = binom3(rvnum);
    for(uint32_t i = 0; i < ltnum; ++i) {
        dst[dst_offset + i] ^= sub[src_offset + i];
    }

    // the rest
    for(uint32_t i = rvnum; i < vnum; ++i) {
        if(!vprod1(v, vnum, i)) {
            continue;
        }

        // reduce to degree 2 monomials
        for(uint32_t j = 1; j < rvnum; ++j) {
            for(uint32_t k = 0; k < j; ++k) {
                dst[midx2(deg, rvnum, k, j)] ^= sub[midx3(deg, vnum, k, j, i)];
            }
        }

        for(uint32_t j = rvnum; j < i; ++j) {
            if(!vprod1(v, vnum, j)) {
                continue;
            }

            // reduce to linear term
            for(uint32_t k = 0; k < rvnum; ++k) {
                dst[midx1(deg, rvnum, k)] ^= sub[midx3(deg, vnum, k, j, i)];
            }

            // reduce to constant term
            for(uint32_t k = rvnum; k < j; ++k) {
                if(vprod1(v, vnum, k)) {
                    dst[midx0(deg, rvnum)] ^= sub[midx3(deg, vnum, k, j, i)];
                }
            }
        }
    }

    // degree 2 monomials
    dst_offset = midx2(deg, rvnum, 0, 1);
    src_offset = midx2(deg, vnum, 0, 1);
    ltnum = binom2(rvnum);
    for(uint32_t i = 0; i < ltnum; ++i) {
        dst[dst_offset + i] ^= sub[src_offset + i];
    }

    // the rest
    for(uint32_t i = rvnum; i < vnum; ++i) {
        if(vprod1(v, vnum, i)) {
            // reduce to linear terms
            for(uint32_t j = 0; j < rvnum; ++j) {
                dst[midx1(deg, rvnum, j)] ^= sub[midx2(deg, vnum, j, i)];
            }

            // reduce to constant
            for(uint32_t j = rvnum; j < i; ++j) {
                if(vprod1(v, vnum, j)) {
                    dst[midx0(deg, rvnum)] ^= sub[midx2(deg, vnum, j, i)];
                }
            }
        }
    }

    // linear terms
    dst_offset = midx1(deg, rvnum, 0);
    src_offset = midx1(deg, vnum, 0);
    for(uint32_t i = 0; i < rvnum; ++i) {
        dst[dst_offset + i] ^= sub[src_offset + i];
    }

    // reduce the rest to constant term
    dst_offset = midx0(deg, rvnum);
    for(uint32_t i = rvnum; i < vnum; ++i) {
        if(vprod1(v, vnum, i)) {
            dst[dst_offset] ^= sub[src_offset + i];
        }
    }

    // constant term
    dst[dst_offset] ^= sub[midx0(deg, vnum)];
}

} /* extern "C" */
