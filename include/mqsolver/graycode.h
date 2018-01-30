/* graycode.h: header file that defines functions for graycode enumeration on
 * GPU
 */

#ifndef __GRAYCODE_H__
#define __GRAYCODE_H__

#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include <cuda_runtime.h>

#include "mq_math.h"
#include "cuda_util.h"
#include "fix.h"


/* ===========================================================================
 *                              GPU architecture
 * ===========================================================================
 */

#define WARP_SIZE               (32) // number of threads in a warp
#define WARP_PER_BLOCK          (4)  // number of warps in a block

/* macro: global_tid
 * usage: compute the global thread id based on blockIdx and threadIdx
 * arguments: void
 * return: global thread id
 */
#define global_tid() \
    (blockIdx.x * blockDim.x + threadIdx.x)

/* macro: warp_tid
 * usage: compute the thread id in a warp based on blockIdx and threadIdx
 * arguments: void
 * return: warp thread id
 */
#define warp_tid() \
    (global_tid() % WARP_SIZE)

/* macro: warp_id
 * usage: compute the global warp id for a thread
 * arguments: void
 * return: the global warp id
 */
#define global_wid() \
    (global_tid() / WARP_SIZE)

/* macro: block_wid
 * usage: compute the warp id in a block
 * arguments: void
 * return: the warp id in a block
 */
#define block_wid() \
    (global_wid() % WARP_PER_BLOCK)

/* storage scheme:
 *      1) the evaluation of f, pf1, pf2, pf3 are stored in global memory
 *      2) last order pf is stored in constant memory
 *      3) use registers to perform Gausssian elimination.
 *
 * Since constant memory holds last order partial derivatives, the size of the
 * last order partial derivs should be less than 64KB.
 */
#define CMEM_SIZE               (64) // KB
#define CMEM_RESV               (4) // reserved for the device, also in KB
#define SLOT_SIZE               (sizeof(uint32_t))
#define MAX_SLOT_NUM            ((CMEM_SIZE - CMEM_RESV) * 1024 / SLOT_SIZE)

extern __device__ __constant__ uint32_t ldderivs[MAX_SLOT_NUM];

/* default max number of solution candidates returned by all GPU threads */
#define MAX_SC_NUM              (0x1U << 15)

/* number of uint32_t in a slot of the mailbox, we need 3 for thread fix vars,
 * kernel fix vars, and linear system solution.
 */
#define MAILBOX_SLOT_INUM       (3)

/* layout of the mailbox:
 *
 *           number of linear systems that don't have enough indepdent eqs.
 *           |
 *           |  total number of solution candidates, including the one that
 *           |  were filtered.
 *           |  |
 *           v  v
 *      [0] [1] [2] [3-5] [6-8] [9-11] ....
 *      
 *      ^            ^      ^     ^
 *      |            |      |     |
 *      |            |      |     3rd slot
 *      |            |      |
 *      |            |      2nd slot
 *      |            |
 *      |            1st slot
 *      |
 *      next slot to write to, initially zero. When the kernel is done, if its
 *      value reaches MAX_FP_NUM, then the threads have produced more candidates
 *      than the mailbox can hold. The correct solution might therefore be
 *      dropped.
 */

/* macro: ldpf3_idx
 * usage: given the index of variables against which partial derivs were taken,
 *      return the index in constant memory
 * arguments: sorted
 *      1) idx1: for x_n, pass 0, for x_{n-1}, pass 1
 *      2) idx2
 *      3) idx3
 * return: the index for the partial derivs in constant memory
 */
#define ldpf3_idx(idx1, idx2, idx3) \
    ((idx1) + binom2(idx2) + binom3(idx3))

/* macro: ldpf4_idx
 * usage: same as ldpf3_idx, for degree 4
 * arguments: sorted
 *      1) idx1: for x_n, pass 0, for x_{n-1}, pass 1
 *      2) idx2
 *      3) idx3
 *      4) idx4
 * return: the index for the partial derivs in constant memory
 */
#define ldpf4_idx(idx1, idx2, idx3, idx4) \
    ((idx1) + binom2(idx2) + binom3(idx3) + binom4(idx4))

/* function: gc_mbsize
 * usage: compute the size of memory block needed for mailbox
 * arguments:
 *      1) max_fpnum: max number of false positives to keep
 * return: the size in bytes
 */
__host__ uint32_t
gc_mbsize(const uint32_t max_fpnum);

/* function: gc_mbcount
 * usage: return the address to the current number of solutions in the mailbox.
 *      The value stored at the address is also the next slot to write to.
 * arguments:
 *      1) mb: poiner to the mailbox
 * return: aforementioned address
 */
__host__ __device__ uint32_t*
gc_mbcount(uint32_t* const mb);

/* function: gc_mbdepc
 * usage: return the address to number of linear system that doesn't have
 *      enough independent eqs to yield a unique solution
 * arguments:
 *      1) mb: poiner to the mailbox
 * return: aforementioned address
 */
__host__ __device__ uint32_t*
gc_mbdepc(uint32_t* const mb);

/* function: gc_mbsnum
 * usage: return the address to number of solution candidates including the
 *      ones that were filtered
 * arguments:
 *      1) mb: poiner to the mailbox
 * return: aforementioned address
 */
__host__ __device__ uint32_t*
gc_mbsnum(uint32_t* const mb);

/* function: gc_mbslot
 * usage: given the index of a slot, compute its address in the mailbox
 * arguments:
 *      1) mb: poiner to the mailbox
 *      2) idx: slot index
 * return: the address to the start of the slot
 */
__host__ __device__ uint32_t*
gc_mbslot(uint32_t* const mb, const uint32_t idx);

/* function: gc_mblsys
 * usage: extract the solution to linear system from one slot of the mailbox
 * arguments:
 *      1) slot: pointer to the slot
 * return: address to the solution to the linear system
 */
__host__ __device__ uint32_t*
gc_mblsys(uint32_t* const mbslot);

/* function: gc_mbtsol
 * usage: extract the value of variables fixed by the thread from one slot
 *      of the mailbox
 * arguments:
 *      1) slot: pointer to the slot
 * return: address to the aforementioned values
 */
__host__ __device__ uint32_t*
gc_mbtsol(uint32_t* const mbslot);

/* function: gc_mbksol
 * usage: extract the value of variables fixed by the kernel from one slot
 *      of the mailbox
 * arguments:
 *      1) slot: pointer to the slot
 * return: address to aforementioned values
 */
__host__ __device__ uint32_t*
gc_mbksol(uint32_t* const mbslot);

/* function: gc_wslot_num
 * usage: compute the number of uint32_t needed for one warp
 * arguments:
 *      1) tfnum: number of variables to fix by a thread
 *      2) knum: number of variables to keep
 *      3) deg: degree of the Macaulay matrix. Can be 3 or 4
 * return: the number of uint32_t needed for a warp
 */
__host__ __device__ uint32_t
gc_wslot_num(const uint32_t tfnum, const uint32_t knum, const uint32_t deg);

/* function: gc_wmemsize
 * usage: compute the size of memory needed for one warp
 * arguments:
 *      1) tfnum: number of variables to fix by a thread
 *      2) knum: number of variables to keep
 *      3) deg: degree of the Macaulay matrix. Can be 3 or 4
 * return: the size of memory in bytes
 */
__host__ uint32_t
gc_wmemsize(const uint32_t tfnum, const uint32_t knum, const uint32_t deg);

/* function: gc_cmemsize
 * usage: compute the size of memory needed for last order partial dervs
 * arguments:
 *      1) tfnum: number of variables to fix by a thread
 *      2) deg: degree of the Macaulay. Can be 3 or 4
 * return: the size of memory needed in bytes
 */
__host__ uint32_t
gc_cmemsize(const uint32_t tfnum, const uint32_t deg);

/* function: gc_smemsize
 * usage: compute the size of shared memory needed
 * arguments:
 *      1) knum: number of variables to keep
 *      2) warp_num: number of warps in a thread block
 * return: the size of memory needed in bytes
 */
__host__ uint32_t
gc_smemsize(const uint32_t knum, const uint32_t warp_num);

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
            const uint32_t kfnum, const uint32_t knum, const uint32_t deg);

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
        const uint32_t kfnum, const uint32_t knum);

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
            bool* const __restrict__ tmp);

#endif  /* __GRAYCODE_H__ */
