/* graycode.cu: implementation of graycode.h
 */

extern "C" {

#include "graycode.h"

/* constant memory for storing last order partial derivatives */
__device__ __constant__ uint32_t ldderivs[MAX_SLOT_NUM];

/* function: gc_mb_size
 * usage: compute the size of memory block needed for mailbox
 * arguments:
 *      1) max_fpnum: max number of false positives to keep
 * return: the size in bytes
 */
__host__ uint32_t
gc_mbsize(const uint32_t max_fpnum) {
    return (3 + MAILBOX_SLOT_INUM * max_fpnum) * sizeof(uint32_t);
}

/* function: gc_mbcount
 * usage: retrieve the current number of solutions in the mailbox. The value is also
 *      the next slot to write to.
 * arguments:
 *      1) mb: poiner to the mailbox
 * return: aforementioned value
 */
__host__ __device__ uint32_t*
gc_mbcount(uint32_t* const mb) {
    return mb;
}

/* function: gc_mbdepc
 * usage: return the address to number of linear system that doesn't have
 *      enough independent eqs to yield a unique solution
 * arguments:
 *      1) mb: poiner to the mailbox
 * return: aforementioned address
 */
__host__ __device__ uint32_t*
gc_mbdepc(uint32_t* const mb) {
    return mb + 1;
}

/* function: gc_mbsnum
 * usage: return the address to number of solution candidates including the
 *      ones that were filtered
 * arguments:
 *      1) mb: poiner to the mailbox
 * return: aforementioned address
 */
__host__ __device__ uint32_t*
gc_mbsnum(uint32_t* const mb) {
    return mb + 2;
}

/* function: gc_mbslot
 * usage: given the index of a slot, compute its address in the mailbox
 * arguments:
 *      1) mb: poiner to the mailbox
 *      2) idx: slot index
 * return: the address to the start of the slot
 */
__host__ __device__ uint32_t*
gc_mbslot(uint32_t* const mb, const uint32_t idx) {
    return mb + 3 + MAILBOX_SLOT_INUM * idx;
}

/* function: gc_mblsys
 * usage: extract the solution to linear system from one slot of the mailbox
 * arguments:
 *      1) slot: pointer to the slot
 * return: address to the solution to the linear system
 */
__host__ __device__ uint32_t*
gc_mblsys(uint32_t* const mbslot) {
    return mbslot;
}

/* function: gc_mbtsol
 * usage: extract the value of variables fixed by the thread from one slot
 *      of the mailbox
 * arguments:
 *      1) slot: pointer to the slot
 * return: address to the aforementioned values
 */
__host__ __device__ uint32_t*
gc_mbtsol(uint32_t* const mbslot) {
    return mbslot + 1;
}

/* function: gc_mbksol
 * usage: extract the value of variables fixed by the kernel from one slot
 *      of the mailbox
 * arguments:
 *      1) slot: pointer to the slot
 * return: address to aforementioned values
 */
__host__ __device__ uint32_t*
gc_mbksol(uint32_t* const mbslot) {
    return mbslot + 2;
}

/* function: gc_boffset
 * usage: compute the number of uint32_t used for one evalution, which is the
 *      offset for one block in the data structures for graycode enumeration
 * arguments:
 *      1) knum: number of variables to keep
 * return: the number of uint32_t used
 */
__device__ static __forceinline__ uint32_t
gc_boffset(const uint32_t knum) {
    return (knum + 1) * WARP_SIZE;
}

/* function: gc_wslot_num
 * usage: compute the number of uint32_t needed for one warp
 * arguments:
 *      1) tfnum: number of variables to fix by a thread
 *      2) knum: number of variables to keep
 *      3) deg: degree of the Macaulay matrix. Can be 3 or 4
 * return: the number of uint32_t needed for a warp
 */
__host__ __device__ uint32_t
gc_wslot_num(const uint32_t tfnum, const uint32_t knum, const uint32_t deg) {
    uint32_t int_num = WARP_SIZE * (knum+1); // f
    int_num += tfnum * WARP_SIZE * (knum+1); // 1st order partial dervs
    int_num += binom2(tfnum) * WARP_SIZE * (knum+1); // 2nd order

    if(4 == deg) {
        int_num += binom3(tfnum) * WARP_SIZE * (knum+1); // 3rd order
    }

    return int_num;
}

/* function: gc_wmemsize
 * usage: compute the size of memory needed for one warp
 * arguments:
 *      1) tfnum: number of variables to fix
 *      2) knum: number of variables to keep
 *      3) deg: degree of the Macaulay matrix. Can be 3 or 4
 * return: the size of memory in bytes
 */
__host__ uint32_t
gc_wmemsize(const uint32_t tfnum, const uint32_t knum, const uint32_t deg) {
    return gc_wslot_num(tfnum, knum, deg) * sizeof(uint32_t);
}

/* function: gc_cmemsize
 * usage: compute the size of memory needed for last order partial dervs
 * arguments:
 *      1) tfnum: number of variables to fix by a thread
 *      2) deg: degree of the Macaulay. Can be 3 or 4
 * return: the size of memory needed in bytes
 */
__host__ uint32_t
gc_cmemsize(const uint32_t tfnum, const uint32_t deg) {
    assert(3 == deg || 4 == deg);
    switch(deg) {
        case 3: return binom3(tfnum) * sizeof(uint32_t);
        case 4: return binom4(tfnum) * sizeof(uint32_t);
        default: return 0;
    }
}

/* function: gc_smemsize
 * usage: compute the size of shared memory needed
 * arguments:
 *      1) knum: number of variables to keep
 *      2) warp_num: number of warps in a thread block
 * return: the size of memory needed in bytes
 */
__host__ uint32_t
gc_smemsize(const uint32_t knum, const uint32_t warp_num) {
    return (knum + 1) * WARP_SIZE * warp_num * sizeof(uint32_t);
}

/* function: gc_cmeminit
 * usage: prepare last order partial derivs on the host
 * arguments:
 *      1) cmem, container for the last order partial derivs
 *      2) sub: the tranposed sub-system with redundant non-linear terms
 *      3) tfnum: number of variables to fix by a thread
 *      4) kfnum: number of variables to fix by a kernel
 *      5) knum: number of variables to keep
 *      6) deg: degree of the partial dervs, can be 3 or 4
 * return: void
 */
__host__ void 
gc_cmeminit(uint32_t* const __restrict__ cmem,
            const uint32_t* const __restrict__ sub, const uint32_t tfnum,
            const uint32_t kfnum, const uint32_t knum, const uint32_t deg) {
    assert(3 == deg || 4 == deg);

    const uint32_t vnum = tfnum + kfnum + knum;

    if(3 == deg) {
        for(uint64_t i = 0; i < tfnum; ++i) {
            for(uint64_t j = i+1; j < tfnum; ++j) {
                for(uint64_t k = j+1; k < tfnum; ++k) {

                    uint32_t idx = midx3(3, vnum, vnum-1-kfnum-k, vnum-1-kfnum-j,
                                         vnum-1-kfnum-i);
                    cmem[ldpf3_idx(i, j, k)] = sub[idx];
                }
            }
        }
        return;
    }

    // degree 4
    for(uint64_t i = 0; i < tfnum; ++i) {
        for(uint64_t j = i+1; j < tfnum; ++j) {
            for(uint64_t k = j+1; k < tfnum; ++k) {
                for(uint64_t l = k+1; l < tfnum; ++l) {

                    uint32_t idx = midx4(vnum, vnum-1-kfnum-l, vnum-1-kfnum-k,
                                         vnum-1-kfnum-j, vnum-1-kfnum-i);
                    cmem[ldpf4_idx(i, j, k, l)] = sub[idx];
                }
            }
        }
    }
}

/* function: gc_warp_data
 * usage: given the starting address of graycode data structures,
 *      return the starting address of the memory block for the warp
 * arguments:
 *      1) gc_data: starting address for graycode enumeration data structures
 *      2) warp_data_size: number of uint32_t used for a warp
 * return: the address of the block for f
 */
__device__ static __forceinline__ uint32_t*
gc_warp_data(uint32_t* const gc_data, const uint32_t warp_data_size) {
    return gc_data + global_wid() * warp_data_size;
}

/* function: gc_warp_f
 * usage: given the starting address of graycode data structures for a warp,
 *      return the starting address of the memory block for f
 * arguments:
 *      1) warp_data: starting address for a warp
 * return: the address of the block for f
 */
__device__ static __forceinline__ uint32_t*
gc_warp_f(uint32_t* const warp_data) {
    return warp_data;
}

/* function: gc_warp_pf1
 * usage: given the starting address of graycode data structures for a warp,
 *      return the starting address of the memory block for a 1st order
 *      partial derivs
 * arguments:
 *      1) wdata: starting address for a warp
 *      2) boffset: number of uint32_t used for one evaluation
 *      3) var_idx: the variable index, should be 0 ~ tfnum-1
 * return: the address of the block for the 1st order partial derivs
 */
__device__ static __forceinline__ uint32_t*
gc_warp_pf1(uint32_t* const wdata, const uint32_t boffset,
            const uint32_t var_idx) {
    return wdata + (1 + var_idx) * boffset;
}

/* function: gc_warp_pf2_start
 * usage: return the starting address of the block for the evaluation of 
 *      2nd order partial derivatives
 * arguments:
 *      1) wdata: starting address for a warp
 *      2) boffset: number of uint32_t used for one evaluation
 *      3) tfnum: number of variables to fix by a thread
 * return: the address
 */
__device__ static __forceinline__ uint32_t*
gc_warp_pf2_start(uint32_t* const wdata, const uint32_t boffset,
                  const uint32_t tfnum) {
    return wdata + (1 + tfnum) * boffset;
}

/* function: gc_warp_pf2
 * usage: return the starting address of evaluation of the 2nd order partial
 *      derivatives against 2 particular variables
 * arguments:
 *      1) wdata_pf2: starting address of 2nd paritial derivs for a warp
 *      2) boffset: number of uint32_t used for one evaluation
 *      3) var1_idx: the smaller variable index, should be 0 ~ tfnum-2
 *      4) var2_idx: the larger variable index, should be 1 ~ tfnum-1
 * return the address as uint32_t*
 */
__device__ static __forceinline__ uint32_t*
gc_warp_pf2(uint32_t* const wdata_pf2, const uint32_t boffset,
            const uint32_t var1_idx, const uint32_t var2_idx) {
    return wdata_pf2 + (binom2(var2_idx) + var1_idx) * boffset;
}

/* function: gc_warp_pf3_start
 * usage: return the starting address of the block for the evaluation of 
 *      3rd order partial derivatives
 * arguments:
 *      1) wdata: starting address for a warp
 *      2) boffset: number of uint32_t used for one evaluation
 *      3) tfnum: number of variables to fix by a thread
 * return: the address
 */
__device__ static __forceinline__ uint32_t*
gc_warp_pf3_start(uint32_t* const wdata, const uint32_t boffset,
                  const uint32_t tfnum) {
    return wdata + (1 + tfnum + binom2(tfnum)) * boffset;
}

/* function: gc_warp_pf3
 * usage: return the starting address of evaluation of the 3rd order partial
 *      derivatives against 3 particular variables
 * arguments:
 *      1) wdata_pf3: starting address of 3rd paritial derivs for a warp
 *      2) boffset: number of uint32_t used for one evaluation, should be
 *              1 if degree of Macaulay matrix is 3
 *      3) var1_idx: the smaller variable index, should be 0 ~ tfnum-3
 *      4) var2_idx: the middle variable index, should be 1 ~ tfnum-2
 *      5) var3_idx: the largest variable index, should be 2 ~ tfnum-1
 * return the address as uint32_t*
 */
__device__ static __forceinline__ uint32_t*
gc_warp_pf3(uint32_t* const wdata_pf3, const uint32_t boffset,
            const uint32_t var1_idx, const uint32_t var2_idx,
            const uint32_t var3_idx) {
    return wdata_pf3 + (binom3(var3_idx) + binom2(var2_idx) + var1_idx) * boffset;
}

/* function: gc_warp_getf
 * usage: return the starting address of evaluation of the target paritial
 *      derivatives against some particular variable
 * arguments:
 *      1) wdata: starting address for a warp
 *      2) wdata_pf2: starting address of 2nd partiall derivs for a warp
 *      3) wdata_pf3: starting address of 3rd paritial derivs for a warp
 *      4) boffset: number of uint32_t used for one evaluation, should be
 *              1 if degree of Macaulay matrix is 3
 *      5) order: order of the paritial derivs
 *      6+) indices: sorted indices for variables, for x_n, pass 0.
 * return: the address as uint32_t*
 */
__device__ static __forceinline__ uint32_t*
gc_warp_getf(uint32_t* const wdata, uint32_t* const wdata_pf2,
             uint32_t* const wdata_pf3, const uint32_t boffset,
             const uint32_t order, const uint32_t idx1,
             const uint32_t idx2, const uint32_t idx3) {
    switch(order) {
        case 0: return gc_warp_f(wdata);
        case 1: return gc_warp_pf1(wdata, boffset, idx1);
        case 2: return gc_warp_pf2(wdata_pf2, boffset, idx1, idx2);
        case 3: return gc_warp_pf3(wdata_pf3, boffset, idx1, idx2, idx3);
        default: return NULL;
    }
}

/* function: gc_warp_smem
 * usage: given the global warp id, compute the starting address of the temporary
 *      storage for the thread in shared memory
 * arguments:
 *      1) smem: address to the start of shared memory
 *      2) knum: number of variables to keep
 * return: aformentioned address
 */
__device__ static __forceinline__ uint32_t*
gc_warp_smem(uint32_t* const smem, const uint32_t knum) {
    return smem + block_wid() * gc_boffset(knum) + warp_tid() * (knum + 1);
}

/* function: gc_init
 * usage: initialize the data structures for graycode enumeration
 * arguments:
 *      1) gc_data: pointer to graycode data structures
 *      2) sub: the tranposed sub-system with redundant non-linear terms
 *      3) deg: degree of the Macaulay matrix
 *      4) vnum: total number of variables in the sub-system
 *      5) tfnum: number of variables to fix by a thread
 *      6) kfnum: number of variables to fix by the kernel. Note that this
 *              means there are 2^kf threads
 *      7) knum: number of variables to keep
 * return: void
 */
__global__ void
gc_init(uint32_t* const __restrict__ gc_data,
        const uint32_t* const __restrict__ sub,
        const uint32_t deg, const uint32_t vnum, const uint32_t tfnum,
        const uint32_t kfnum, const uint32_t knum) {
    assert(3 == deg || 4 == deg);
    assert(vnum == (tfnum + kfnum + knum));

    const uint32_t kv = global_tid();
    uint32_t* const wdata = gc_warp_data(gc_data, gc_wslot_num(tfnum, knum, deg));
    const uint32_t boffset = gc_boffset(knum);
    uint32_t* const wdata_pf2 = gc_warp_pf2_start(wdata, boffset, tfnum);
    uint32_t* const wdata_pf3 = gc_warp_pf3_start(wdata, boffset, tfnum);

    // evaluate the sub-system by setting x_{n-kf+1} ~ x_n according to global
    // thread id and the rest tf variables x_{k+1} ~ x_{n-kf} to zero.
    fix_f(gc_warp_f(wdata), sub, deg, vnum, knum, kfnum, global_tid());

    // evaluate paritial derivs of the sub-system
    if(3 == deg) {
        
        for(uint32_t i = 0; i < tfnum; ++i) {
            uint32_t input = 0x0U;
            if(0 != i) {
                input ^= 0x1U << (i-1);
            }

            // init f \par x_? with the input
            // NOTE: i represents x_{var_num-i}, which has index var_num-1-i
            fix_pf1(gc_warp_pf1(wdata, boffset, i), sub, 3, vnum, knum,
                    kfnum, input, kv, vnum-1-kfnum-i);

            for(uint32_t j = i+1; j < tfnum; ++j) {
                uint32_t input2 = input;
                if( i+1 != j) {
                    input2 ^= 0x1U << (j-1);
                }

                // init f \par x_? x_? with the input
                fix_pf2(gc_warp_pf2(wdata_pf2, boffset, i, j), sub, 3, vnum,
                        knum, kfnum, input2, kv, vnum-1-kfnum-j, vnum-1-kfnum-i);
                
                // degree 3 partial derivs are stored in constant memory
            }
        }

        return;
    }

    // degree 4
    for(uint32_t i = 0; i < tfnum; ++i) {
        uint32_t input = 0x0U;
        if(0 != i) {
            input ^= 0x1U << (i-1);
        }

        // init f \par x_? with the input
        // NOTE: i represents x_{var_num-i}, which has index var_num-1-i
        fix_pf1(gc_warp_pf1(wdata, boffset, i), sub, 4, vnum, knum,
                kfnum, input, kv, vnum-1-kfnum-i);

        for(uint32_t j = i+1; j < tfnum; ++j) {
            uint32_t input2 = input;
            if( i+1 != j) {
                input2 ^= 0x1U << (j-1);
            }

            // init f \par x_? x_? with the input
            fix_pf2(gc_warp_pf2(wdata_pf2, boffset, i, j), sub, 4, vnum,
                    knum, kfnum, input2, kv, vnum-1-kfnum-j, vnum-1-kfnum-i);
            
            for(uint32_t k = j+1; k < tfnum; ++k) {
                uint32_t input3 = input2;
                if( j+1 != k) {
                    input3 ^= 0x1U << (k-1);
                }

                // init f \par x_? x_? x_? with the input
                fix_pf3(gc_warp_pf3(wdata_pf3, boffset, i, j, k), sub, vnum,
                        knum, kfnum, input3, kv, vnum-1-kfnum-k, vnum-1-kfnum-j,
                        vnum-1-kfnum-i);

                // degree 4 partial derivs are stored in constant memory
            }
        }
    }
}

/* count the number of trailing zeros, the input c must not be zero */
__device__ static __forceinline__ uint32_t
ctz(uint32_t c) {
    //return __ffs(c) - 1;
    return __clz(__brev(c));
}

/* function: gc_update
 * usage: update the data structures for graycode enumeration
 * arguments:
 *      1) wdata: starting address for a warp
 *      2) wdata_pf2: starting address for pf2 in a warp
 *      3) wdata_pf3: starting address for pf3 in a warp
 *      4) deg: degree of the Macaulay matrix
 *      5) knum: number of variables to keep
 *      6) warp_tid: thread id in a warp
 *      7) c: counter for the enumeration, should not be 0
 * return: void
 */
__device__ static __forceinline__ void
gc_update(uint32_t* const __restrict__ wdata,
          uint32_t* const __restrict__ wdata_pf2,
          uint32_t* const __restrict__ wdata_pf3,
          const uint32_t deg, const uint32_t knum,
          const uint32_t warp_tid, const uint32_t c) {
    assert(0 != c);

    uint32_t idx[4]; // at most degree 4

    const uint32_t pcount = __popc(c);
    const uint32_t boffset = gc_boffset(knum);
    uint32_t bound = (pcount > deg) ? deg : pcount;

    idx[0] = ctz(c);
    //uint32_t rmask = ~0x0U ^ (0x1U << idx[0]);
    for(uint32_t i = 1; i < bound; ++i) {
        //idx[i] = ctz(c & rmask);
        //rmask ^= 0x1U << idx[i];
        idx[i] = ctz(c >> (idx[i-1]+1)) + idx[i-1]+1;
    }

    if(deg == bound) { // update with last partial derivs

        uint32_t* const dest = (3 == deg) ?
                               gc_warp_pf2(wdata_pf2, boffset, idx[0], idx[1]) :
                               gc_warp_pf3(wdata_pf3, boffset, idx[0], idx[1], idx[2]);

        const uint32_t src_idx = (3 == deg) ?
                                 ldpf3_idx(idx[0], idx[1], idx[2]) :
                                 ldpf4_idx(idx[0], idx[1], idx[2], idx[3]);

        eterm_lsys(dest, warp_tid, knum) ^= ldderivs[src_idx];
        --bound;
    }

    for(uint32_t i = bound; i >= 1; --i) {
        // when the number of variables to fix is less than 2,
        // there is no 2nd (and/or 1st) order partial derivs
        uint32_t* const dest = gc_warp_getf(wdata, wdata_pf2, wdata_pf3, boffset,
                                           i-1, idx[0], idx[1], idx[2]);
        uint32_t* const src = gc_warp_getf(wdata, wdata_pf2, wdata_pf3, boffset,
                                           i, idx[0], idx[1], idx[2]);
        // xor k+1 terms
        for(uint32_t j = 0; j < knum+1; ++j) {
            eterm_lsys(dest, warp_tid, j) ^= eterm_lsys(src, warp_tid, j);
        }
    }
}

/* function: gc_copy_lsys
 * usage: copy the linear system for a thread from the interleaved global
 *      memory to a consecutive memory block
 * arguments:
 *      1) mem: address of the memory block
 *      2) f: address to the interleaved global memory where lsys is stored
 *      3) knum: number of variables to keep
 *      4) warp_tid: thread id in a warp
 * return: void
 */
__device__ static __forceinline__ void
gc_copy_lsys(uint32_t* const mem, const uint32_t* const f,
             const uint32_t knum, const uint32_t warp_tid) {
    for(uint32_t i = 0; i < knum+1; ++i) {
        mem[i] = eterm_lsys(f, warp_tid, i);
    }
}

/* function: gc_check_lsys
 * usage: given a transposed linear system, check if it is solvable or not.
 *      The system is stored as an array of uint32_t where each uint32_t holds
 *      the same column for all rows.
 * arguments:
 *      1) sys: the transposed linear system
 *      2) col_num: number of columns
 * return: true if solvable, otherwise false
 */
__device__ static __forceinline__ bool
gc_check_lsys(uint32_t* const __restrict__ sys, const uint32_t col_num) {
    uint32_t rmask = ~0x0U;

    for(uint32_t i = 0; i < col_num-1; ++i) {
        const uint32_t tmp = sys[i] & rmask;
        const uint32_t sf = (!tmp) ? 0x0U : ~0x0U;
        const uint32_t piv = ctz(tmp);
        rmask ^= (0x1U << piv) & sf;
        const uint32_t mask = ( sys[i] ^ (0x1U << piv) ) & sf;

        // for computing the solution after checking
        sys[i] ^= mask;

        for(uint32_t j = i+1; j < col_num; ++j) {
            uint32_t b = ( (sys[j] >> piv) & 0x1U ) ? ~0x0U : 0x0U;
            sys[j] ^= mask & b;
        }
    }

    return !(sys[col_num-1] & rmask);
}

/* function: gc_gjelim
 * usage: given a transposed linear system, reduce it to reduced row echelon
 *      form using Gauss-Jordan elimination. The system is stored as an array
 *      of uint32_t where each uint32_t holds the same column for all rows.
 * arguments:
 *      1) sys: the transposed linear system
 *      2) col_num: number of columns
 * return: void
 */
__device__ static __forceinline__ void
gc_gjelim(uint32_t* const __restrict__ sys, const uint32_t col_num) {
    for(uint32_t i = 0; i < col_num; ++i) {
        
        // NOTE: one can use this flag to avoid branching
        //bool sf = sys[i] >> i; // skip flag
        //uint32_t skip_mask = (sf ^ 0x1U) - 0x1U; // sign extend sf
        // all 1's if don't skip, 0x0UL if skip

        if( !(sys[i] >> i) ) {
            continue; // singular
        }

        // find pivot position, which is the number of trailing zeros
        const uint32_t piv = ctz(sys[i] >> i) + i;

        //const uint32_t mask = (sys[i] ^ (0x1UL << piv)) & skip_mask;
        const uint32_t mask = sys[i] ^ (0x1U << piv);
        sys[i] ^= mask;

        // swap rows, i.e. swap all the ith bits with piv-th bits
        for(uint32_t j = 0; j < i+1; ++j) {
            // i-th bit ^ piv-th bit
            uint32_t tmp = ((sys[j] >> i) ^ (sys[j] >> piv)) & 0x1U;
            //sys[j] ^= ((tmp << i) | (tmp << piv)) & skip_mask;
            sys[j] ^= (tmp << i) | (tmp << piv);
        }

        // row reduction
        for(uint32_t j = i+1; j < col_num; ++j) {
            uint32_t b = (sys[j] >> piv) & 0x1U;
            //sys[j] ^= mask & (~b + 1);
            sys[j] ^= mask & ((b ^ 0x1U) - 0x1U); // sign extend b

            // swap rows, i.e. swap all the ith bits with piv-th bits
            // i-th bit ^ piv-th bit
            uint32_t tmp = ((sys[j] >> i) ^ (sys[j] >> piv)) & 0x1U;
            //sys[j] ^= ((tmp << i) | (tmp << piv)) & skip_mask;
            sys[j] ^= (tmp << i) | (tmp << piv);
        }
    }
}

/* function: gc_solve_lsys
 * usage: given a transposed linear system, perform Gaussian-Jordan elimination
 *      to solve the linear system. The system is stored as an array of
 *      uint32_t where each uint32_t holds the column for all rows. Note that
 *      the system is modified in place.
 * arguments:
 *      1) sys: the transposed linear system
 *      2) col_num: number of columns
 * return: true if the system has a unique or multiple solutions, false
 *      otherwise
 */
__device__ static __forceinline__ bool
gc_solve_lsys(uint32_t* const __restrict__ sys, const uint32_t col_num) {
    gc_gjelim(sys, col_num);
    return !((sys[col_num-1] >> (col_num-1)) & 0x1UL);
}

/* function: gc_extract_sol
 * usage: extract solution of a linear system from a system that has been
 *      processed by gc_check_lsys. Must have no missing pivots, i.e. no
 *      column of the system is all zero.
 * arguments:
 *      1) sys: processed linear system
 *      2) vnum: number of variables in the linear system, should be one
 *              less than the number of columns
 * return: the solution
 */
__device__ static __forceinline__ uint32_t
gc_extract_sol(const uint32_t* const __restrict__ sys, const uint32_t vnum) {
    uint32_t sol = 0x0U;

    for(uint32_t i = 0; i < vnum; ++i) {
        sol |= ( (sys[vnum] >> ctz(sys[i])) & 0x1U ) << i;
    }

    return sol;
}

/* function: gc_check_dep
 * usage: check a linear system that has been processed by gc_check_lsys to
 *      see if there were linear dependent eqs.
 * arguments:
 *      1) sys: processed linear system
 *      2) vnum: number of variables in the linear system, should be one
 *              less than the number of columns
 * return: 1 if there were, otherwise 0
 */
__device__ static __forceinline__ uint32_t
gc_check_dep(const uint32_t* const __restrict__ sys, const uint32_t vnum) {

    for(uint32_t i = 0; i < vnum; ++i) {
        if(!sys[i]) { // all zero
            return 0x1U;
        }
    }

    return 0x0U;
}

/* subroutine of gc_bfsearch, verify a solution candidate with 32 filter eqs */
__device__ static __forceinline__ bool
gc_filter_sol(const uint32_t* const __restrict__ mq_eqs, const uint32_t lsys_sol,
              const uint32_t count, const uint32_t kf_sol, const uint32_t sfval,
              const uint32_t knum, const uint32_t tfnum, const uint32_t kfnum,
              const uint32_t sfnum, bool* const __restrict__ tmp) {

    const uint32_t gc = count ^ (count >> 1);
    const uint32_t front_num = knum + tfnum;
    const uint32_t svnum = front_num + kfnum;
    const uint32_t vnum = svnum + sfnum;

    // TODO: use registers?
    bool* sol = tmp + kf_sol * vnum;

    for(uint32_t i = 0; i < knum; ++i) {
        sol[i] = (lsys_sol >> i) & 0x1U;
    }
    
    for(uint32_t i = knum; i < front_num; ++i) {
        sol[i] = (gc >> (front_num-1-i)) & 0x1U;
    }

    for(uint32_t i = front_num; i < svnum; ++i) {
        sol[i] = (kf_sol >> (svnum-1-i)) & 0x1U;
    }

    for(uint32_t i = svnum; i < vnum; ++i) {
        sol[i] = (sfval >> (vnum-1-i)) & 0x1U;
    }

    uint64_t pace = 0;
    uint32_t res = 0x0U;

    // degree 2 terms
    for(uint32_t i = 1; i < vnum; ++i) {
        if(!sol[i]) {
            pace += i;
            continue;
        }

        for(uint32_t j = 0; j < i; ++j) {
            if(sol[j]) {
                res ^= mq_eqs[pace];
            }

            ++pace;
        }
    }

    // linear terms
    for(uint32_t i = 0; i < vnum; ++i) {
        if(sol[i]) {
            res ^= mq_eqs[pace];
        }
        ++pace;
    }

    // constant term
    return !(res ^ mq_eqs[pace]);
}

/* function: gc_bfsearch
 * usage: split the search space evenly into subspaces by fixing variables at
 *      kernel level and use each thread to bruteforce the subspace.
 * arguments:
 *      1) gc_data: pointer to graycode data structures
 *      2) deg: degree of the Macaulay matrix
 *      3) tfnum: number of variables to fix by a thread
 *      4) kfnum: number of variables to fix by a kernel
 *      5) sfnum: number of variables fixed in the sub-system
 *      6) sfval: store the values for variables fixed in the sub-system
 *      7) mb: pointer to the mailbox
 *      8) mbsize: max number of solutions for the mailbox
 *      9) mq_eqs: 32 eqs from the MQ system stored in column-wise format
 *      10) tmp: temporary storage for verifying a solution candidate
 * return: void
 */
__global__ void
gc_bfsearch(uint32_t* const __restrict__ gc_data, const uint32_t deg,
            const uint32_t tfnum, const uint32_t kfnum, const uint32_t sfnum,
            const uint32_t sfval, uint32_t* const __restrict__ mb,
            const uint32_t mbsize, const uint32_t* const __restrict__ mq_eqs,
            bool* const __restrict__ tmp) {
    assert(tfnum < 32);

    // temporary storage/registers for performing Gaussian elimination
    #include "gc_decl_lsys.def"

    const uint32_t boffset = gc_boffset(KNUM);
    const uint32_t warp_tid = warp_tid();
    uint32_t* const wdata = gc_warp_data(gc_data, gc_wslot_num(tfnum, KNUM, deg));
    uint32_t* const wdata_pf2 = gc_warp_pf2_start(wdata, boffset, tfnum);
    uint32_t* const wdata_pf3 = gc_warp_pf3_start(wdata, boffset, tfnum);

    uint32_t count = 0;
    const uint32_t sub_size = 0x1U << tfnum;

    while(count < sub_size) {
        // copy the linear system into registers because Gaussian-Jordan
        // elimination modifies the data in place.
        #include "gc_copy_lsys.def"

        bool solvable;
        #include "gc_check_lsys.def"

        if(solvable) { // found a potential solution

#ifdef GC_DEPC_LSYS

            uint32_t dep = 0x0U;
            #include "gc_dep_lsys.def"

            if(dep) { // there are dependent eqs in sub-sys
                atomicAdd(gc_mbdepc(mb), 1);
                // TODO: dump the linear system for further checking
            }

#endif

            atomicAdd(gc_mbsnum(mb), 1);
            uint32_t sol = 0x0U;
            #include "gc_extract_sol.def"

            // apply the filter
            if(gc_filter_sol(mq_eqs, sol, count, global_tid(), sfval,
                             KNUM, tfnum, kfnum, sfnum, tmp)) {
                const uint32_t slot_idx = atomicAdd(gc_mbcount(mb), 1);
                if(mbsize > slot_idx) {
                    uint32_t* const slot = gc_mbslot(mb, slot_idx);
                    *gc_mblsys(slot) = sol;
                    *gc_mbtsol(slot) = count;
                    *gc_mbksol(slot) = global_tid();
                }
            }
        }

        gc_update(wdata, wdata_pf2, wdata_pf3, deg, KNUM, warp_tid, ++count);
    }
}

} /* extern "C" */
