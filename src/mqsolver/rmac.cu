/* rmac.c: implementation of structure RMac
 */

extern "C" {

#include "rmac.h"

/* function: drm_slots_num
 * usage: compute how many slots are needed
 * argument: mnum, number of monomials
 * return: the number of slots
 */
__host__ static __forceinline__  uint64_t
drm_slots_num(const uint64_t mnum) {
    return (1 + ((mnum - 1) / 64));
}

/* function: drm_slot_idx
 * usage: given the monomial index, compute the index of the slot where it locates
 * arguments:
 *      1) col_idx: index of the monomial
 * return: the slot index
 */
__device__ __host__ static __forceinline__ uint32_t
drm_slot_idx(const uint32_t col_idx) {
    return col_idx >> 6;
}

/* function: drm_slot_offset
 * usage: given the monomial index, compute the offset from RHS into slot where
 *      it belongs
 * arguments:
 *      1) col_idx: index of the monomial
 * return: the slot index
 */
__device__ __host__ static __forceinline__ uint32_t
drm_slot_offset(const uint32_t col_idx) {
    return 63 - (col_idx & 0x3FUL); // modulo 64, then reverse
}

/* function: drm_at
 * usage: return the monomial at the given index
 * arguments:
 *      1) row: a pointer to the dense row
 *      2) col_idx: index of the monomial, 0 is first monomial is grlex
 * return: void
 */
__device__ __host__ static __forceinline__ uint32_t
drm_at(uint64_t* row, const uint32_t col_idx) {
    return (row[drm_slot_idx(col_idx)] >> drm_slot_offset(col_idx)) & 0x1UL;
}

/* function: rmac_create
 * usage: create a container for sub-system
 * arguments:
 *      1) eq_num: number of eqs
 *      2) term_num: number of terms
 * return: pointer to the newly created structure RMac
 */
__host__ RMac*
rmac_create(const uint64_t eq_num, const uint64_t term_num) {
    RMac* s = NULL;
    uint64_t slot_num = drm_slots_num(term_num);

    s = (RMac*) malloc(sizeof(RMac) + sizeof(uint64_t*) * eq_num);
    if(!s) {
        EXIT_WITH_MSG("[!] failed to allocate memory for RMac\n");
    }

    s->eq_num = eq_num;
    s->term_num = term_num;
    s->slot_num = slot_num;
    s->mem = (uint64_t*) calloc(eq_num * slot_num, sizeof(uint64_t));
    if(!s->mem) {
        EXIT_WITH_MSG("[!] failed to allocate memory for RMac\n");
    }
    for(uint64_t i = 0; i < eq_num; ++i) {
        s->eqs[i] = s->mem + slot_num * i;
    }
    
    return s;
}

/* function: rmac_free
 * usage: free structure RMac
 * arguments:
 *      1) sys: pointer to structure RMac
 * return: void
 */
__host__ void
rmac_free(RMac* sys) {
    if(sys->mem) {
        free(sys->mem);
    }
    free(sys);
}

/* function: rmac_memsize
 * usage: given the number of eqs and terms, compute the memory needed for rmac
 * arguments:
 *      1) eq_num: number of eqs
 *      2) term_num: number of terms
 * return: size in bytes
 */
__host__ uint64_t
rmac_memsize(const uint64_t eq_num, const uint64_t term_num) {
    uint64_t memsz = sizeof(RMac) + eq_num * sizeof(uint64_t*);
    return memsz + eq_num * drm_slots_num(term_num) * sizeof(uint64_t);
}

/* function: rmac_row
 * usage: given the row index, return the starting address for the row in RMac
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) idx: row index
 * return: the address
 */
__host__ uint64_t*
rmac_row(RMac* const sys, const uint64_t idx) {
    return sys->eqs[idx];
}

/* function: rmac_at
 * usage: given the row and column index, return the entry
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) ridx: row index
 *      3) cidx: column index
 * return: the entry
 */
__host__ bool
rmac_at(RMac* const sys, const uint64_t ridx, const uint64_t cidx) {
    return drm_at(rmac_row(sys, ridx), cidx);
}

/* function: rmac_zrow
 * usage: check if a row is all zeros
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) row_idx: row index
 * return: true if the row is all zero, otherwise false
 */
__host__ bool
rmac_zrow(RMac* sys, const uint64_t row_idx) {

    uint64_t* row = rmac_row(sys, row_idx);
    for(uint64_t i = 0; i < sys->slot_num; ++i) {
        if( row[i] ) {
            return false;
        }
    }

    return true;
}

/* function: rmac_perm_rows
 * usage: permute the rows according to an index map
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) idx_map: row index map
 * return: void
 */
__host__ void
rmac_perm_rows(RMac* sys, const uint32_t* const idx_map) {
    for(uint64_t i = 0; i < sys->eq_num; ++i) {
        sys->eqs[i] = sys->mem + idx_map[i] * sys->slot_num;
    }
}

/* subroutine of rmac_elim. Initialize data structure for rmac_elim */
static __global__ void
rmac_elim_init(uint64_t* const sys, uint64_t** rows, uint32_t* row_indices,
               const uint32_t eq_num, uint32_t slot_num) {

    const uint32_t tid = global_tid();
    if(tid < eq_num) {
        rows[tid] = sys + tid * slot_num;
        row_indices[tid] = tid;
    }
}

/* subroutine of rmac_elim. Find all the rows below the target pivot row
 * whose pivot monomial is 1. The uint32_t that rcount points to must be
 * initialized to zero.
 */
static __global__ void
find_pvt_rows(uint64_t** __restrict__ rows, const uint32_t eq_num,
              const uint32_t start, const uint32_t col_idx,
              uint32_t* const __restrict__ indices,
              uint32_t* const __restrict__ rcount) {

    const uint32_t tid = global_tid();
    if(start <= tid && tid < eq_num) {
        uint64_t* row = rows[tid];
            
        // check the target monomial in the row
        if( drm_at(row, col_idx) ) {
            indices[atomicAdd(rcount, 1)] = tid;
        }
    }
}

/* subroutine of rmac_elim. Swap two rows and update the row indices
 * accordingly.
 */
static __global__ void
swap_rows(uint64_t** __restrict__ rows, const uint32_t idx1,
          const uint32_t* const __restrict__ idx2_addr,
          uint32_t* const __restrict__ row_indices) {

    if(global_tid() == 0) {
        uint32_t idx2 = *idx2_addr;
        uint64_t* tmp = rows[idx1];
        rows[idx1] = rows[idx2];
        rows[idx2] = tmp;

        uint32_t idx = row_indices[idx1];
        row_indices[idx1] = row_indices[idx2];
        row_indices[idx2] = idx;
    }
}

/* subroutine of rmac_elim. Perform row reduction on all non-pivot rows */
static __global__ void
reduc_rows(uint64_t** __restrict__ rows, const uint32_t slot_num,
           const uint32_t start, uint32_t* const __restrict__ indices,
           const uint32_t rcount) {

    extern __shared__ uint64_t smem[];
    const uint32_t tid = global_tid();

    if(start + tid < slot_num) {
        uint64_t* src = rows[indices[0]];

        // copy src into shared memory
        smem[threadIdx.x] = src[start + tid];

        // no need to sync

        // for each rows below
        for(uint32_t i = 1; i < rcount; ++i) {
            uint64_t* dst = rows[indices[i]];
            dst[start + tid] ^= smem[threadIdx.x];
        }
    }
}

/* function: rmac_elim_gpu
 * usage: perform Gaussian elimination on the reduced Macaulay matrix with GPU
 * arguments:
 *      1) sys: pointer to the memory block that store the reduced Macaulay
 *              matrix
 *      2) rows: points to the starting addresses of rows in the reduced Macaulay
 *              matrix
 *      3) row_indices: used to track the row indices after swapping
 *      4) eq_num: numner of eqations in the reduced Macaulay matrix
 *      5) slot_num: number of uint32_t used for each row
 *      6) col_num: number of columns to perform Gaussian elimination
 *      7) indices: container for the indices of pivot row candidates
 *      8) rcount: point to device memory used to stored the number of pivot row
 *              candidates
 *      9) host_rcount: same as above except it's on the host
 * return: void
 */
__host__ void
rmac_elim_gpu(uint64_t* const sys, uint64_t** rows, uint32_t* row_indices,
              const uint32_t eq_num, const uint32_t slot_num,
              const uint32_t col_num, uint32_t* const __restrict__ indices,
              uint32_t* const __restrict__ rcount,
              uint32_t* const __restrict__ host_rcount) {

    uint32_t warp_num = 1 + (eq_num - 1) / WARP_SIZE;
    uint32_t blk_size = 4;
    uint32_t blk_num = 1 + (warp_num - 1) / blk_size;
    blk_size *= WARP_SIZE;

    CUDA_CHECK(cudaFuncSetCacheConfig(reduc_rows, cudaFuncCachePreferShared));
    CUDA_CHECK(cudaFuncSetSharedMemConfig(reduc_rows,
                                          cudaSharedMemBankSizeEightByte));

    rmac_elim_init<<<blk_num, blk_size>>>(sys, rows, row_indices, eq_num, slot_num);
    CUDA_CHECK(cudaPeekAtLastError());

    for(uint32_t i = 0; i < col_num; ++i) {
        warp_num = 1 + (eq_num - 1) / WARP_SIZE;
        blk_size = 32;
        blk_num = 1 + (warp_num - 1) / blk_size;
        blk_size *= WARP_SIZE;

        CUDA_CHECK(cudaMemset(rcount, 0x0, sizeof(uint32_t)));
        find_pvt_rows<<<blk_num, blk_size>>>(rows, eq_num, i, i, indices, rcount);
        CUDA_CHECK(cudaPeekAtLastError());
        CUDA_CHECK(cudaMemcpy(host_rcount, rcount, sizeof(uint32_t),
                              cudaMemcpyDeviceToHost));
  
        if( ! *host_rcount ) { // singular
            continue;
        }

        // adjust launch config according to rcount and terms to xor
        uint32_t start = drm_slot_idx(i);
        uint32_t rsnum = slot_num - start;
        warp_num = 1 + (rsnum - 1) / WARP_SIZE;
        blk_size = 32;
        blk_num = 1 + (warp_num - 1) / blk_size;
        blk_size *= WARP_SIZE;
        uint32_t smem_size = sizeof(uint64_t) * blk_size;

        reduc_rows<<<blk_num, blk_size, smem_size>>>(rows, slot_num, start, indices,
                                                     *host_rcount);
        CUDA_CHECK(cudaPeekAtLastError());

        swap_rows<<<1, 1>>>(rows, i, indices, row_indices);
        CUDA_CHECK(cudaPeekAtLastError());
    }
}

/* function: rmac_xor_after
 * usage: perform xor on two rows after a certain bit index and stored the
 *      result back to the first row
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) r1: starting address of the first row
 *      3) r2: starting address of the second row
 *      4) idx: bit index
 * return: void
 */
__host__ static inline void
rmac_xor_after(RMac* const sys, uint64_t* const __restrict__ r1,
               const uint64_t* const __restrict__ r2,
               const uint64_t idx) {
    uint64_t start = drm_slot_idx(idx);
    uint64_t end = sys->slot_num;
    uint64_t round = (end - start) / 4;
    uint64_t bound = start + round * 4;

#ifdef MAC_DROW_SSE

    for(uint64_t i = start; i < bound; i += 4) {
        __m128i src1 = _mm_loadu_si128((__m128i*) &r2[i+0]);
        __m128i dst1 = _mm_loadu_si128((__m128i*) &r1[i+0]);
        __m128i src2 = _mm_loadu_si128((__m128i*) &r2[i+2]);
        __m128i dst2 = _mm_loadu_si128((__m128i*) &r1[i+2]);

        _mm_storeu_si128((__m128i*)&r1[i+0], _mm_xor_si128(dst1, src1));
        _mm_storeu_si128((__m128i*)&r1[i+2], _mm_xor_si128(dst2, src2));
    }

#else

    for(uint64_t i = start; i < bound; i += 4) {
        r1[i] ^= r2[i];
        r1[i+1] ^= r2[i+1];
        r1[i+2] ^= r2[i+2];
        r1[i+3] ^= r2[i+3];
    }

#endif

    for(uint64_t i = bound; i < end; ++i) {
        r1[i] ^= r2[i];
    }
}

/* arguments for worker thread of rmac_elim_cpu */
typedef struct {
    RMac* sys;
    uint32_t* indices;
    uint32_t* local_indices;
    pthread_mutex_t* index_lock;
    uint64_t* index_offset;
    uint64_t start;
    uint64_t end;
    uint64_t col_idx;
} fpvt_arg;

/* worker thread of rmac_elim_cpu */
static void
find_pvt_rows_cpu(void* dummy) {
    fpvt_arg* arg = (fpvt_arg*) dummy;

    uint64_t count = 0;
    for(uint64_t i = arg->start; i < arg->end; ++i) {
        uint64_t* const row = rmac_row(arg->sys, i);
        if( drm_at(row, arg->col_idx) ) {
            arg->local_indices[count++] = i;
        }
    }

    if(count) {
        // copy the local indices into global one
        pthread_mutex_lock(arg->index_lock);

        uint64_t offset = *(arg->index_offset);
        *(arg->index_offset) += count;

        pthread_mutex_unlock(arg->index_lock);

        memcpy(arg->indices + offset, arg->local_indices, sizeof(uint32_t) * count);
    }
}

/* arguments for worker thread of rmac_elim_cpu */
typedef struct {
    RMac* sys;
    const uint32_t* indices;
    uint64_t start;
    uint64_t end;
    uint64_t col_idx;
} reduc_arg;

/* worker thread of rmac_elim_cpu */
static void
reduce_rows_cpu(void* dummy) {
    reduc_arg* arg = (reduc_arg*) dummy;

    const uint64_t* const pvt = rmac_row(arg->sys, arg->indices[0]);
    for(uint64_t i = arg->start; i < arg->end; ++i) {
        uint64_t* const dst = rmac_row(arg->sys, arg->indices[i]);
        rmac_xor_after(arg->sys, dst, pvt, arg->col_idx);
    }
}

/* function: rmac_elim_cpu
 * usage: perform Gaussian elimination on the sub-system
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) col_num: number of columns to perform elimination
 *      3) tpool: pointer to structure threadpool_t
 * return: void
 */
__host__ void
rmac_elim_cpu(RMac* const sys, const uint64_t col_num,
              threadpool_t* const tpool) {

    uint64_t index_offset = 0;
    uint32_t* indices = SMALLOC(uint32_t, sys->eq_num);
    fpvt_arg* fpvt_argpool = SMALLOC(fpvt_arg, tpool->thread_count);
    pthread_mutex_t index_lock = PTHREAD_MUTEX_INITIALIZER;
    reduc_arg* rdc_argpool = SMALLOC(reduc_arg, tpool->thread_count);
    uint32_t* local_indices = SMALLOC(uint32_t, sys->eq_num * tpool->thread_count);

    for(uint64_t i = 0; i < col_num; ++i) {
        index_offset = 0;
        uint64_t chunk = (sys->eq_num - i) / tpool->thread_count;

        // find pivot row candidates
        for(int tid = 0; tid < tpool->thread_count-1; ++tid) {
            fpvt_argpool[tid].sys = sys;
            fpvt_argpool[tid].indices = indices;
            fpvt_argpool[tid].local_indices = local_indices + tid * sys->eq_num;
            fpvt_argpool[tid].index_offset = &index_offset;
            fpvt_argpool[tid].index_lock = &index_lock;
            fpvt_argpool[tid].start = i + tid * chunk;
            fpvt_argpool[tid].end = i + (tid + 1) * chunk;
            fpvt_argpool[tid].col_idx = i;

            if(threadpool_add(tpool, find_pvt_rows_cpu, fpvt_argpool + tid, 0)) {
                EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
            }
        }

        fpvt_argpool[tpool->thread_count-1].sys = sys;
        fpvt_argpool[tpool->thread_count-1].indices = indices;
        fpvt_argpool[tpool->thread_count-1].local_indices = local_indices +
                                            (tpool->thread_count-1) * sys->eq_num;
        fpvt_argpool[tpool->thread_count-1].index_offset = &index_offset;
        fpvt_argpool[tpool->thread_count-1].index_lock = &index_lock;
        fpvt_argpool[tpool->thread_count-1].start = i + (tpool->thread_count-1) * chunk;
        fpvt_argpool[tpool->thread_count-1].end = sys->eq_num;
        fpvt_argpool[tpool->thread_count-1].col_idx = i;

        if(threadpool_add(tpool, find_pvt_rows_cpu,
                          fpvt_argpool + tpool->thread_count-1, 0)) {
            EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
        }
        
        if(threadpool_join(tpool, 0)) {
            EXIT_WITH_MSG("[!] failed to join thread pool\n");
        }

        if(!index_offset) { // singular
            continue;
        }

        // perform row reduction
        chunk = (index_offset-1) / tpool->thread_count;
        for(int tid = 0; tid < tpool->thread_count-1; ++tid) {
            rdc_argpool[tid].sys = sys;
            rdc_argpool[tid].indices = indices;
            rdc_argpool[tid].start = 1 + tid * chunk;
            rdc_argpool[tid].end = 1 + (tid + 1) * chunk;
            rdc_argpool[tid].col_idx = i;

            if(threadpool_add(tpool, reduce_rows_cpu, rdc_argpool + tid, 0)) {
                EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
            }
        }

        rdc_argpool[tpool->thread_count-1].sys = sys;
        rdc_argpool[tpool->thread_count-1].indices = indices;
        rdc_argpool[tpool->thread_count-1].start = 1 + (tpool->thread_count-1) * chunk;
        rdc_argpool[tpool->thread_count-1].end = index_offset;
        rdc_argpool[tpool->thread_count-1].col_idx = i;

        if(threadpool_add(tpool, reduce_rows_cpu,
                          rdc_argpool + tpool->thread_count-1, 0)) {
            EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
        }

        if(threadpool_join(tpool, 0)) {
            EXIT_WITH_MSG("[!] failed to join thread pool\n");
        }

        // swap pivot row into final position
        uint64_t* tmp = sys->eqs[i];
        sys->eqs[i] = sys->eqs[indices[0]];
        sys->eqs[indices[0]] = tmp;
    }

    SFREE(local_indices);
    SFREE(fpvt_argpool);
    SFREE(indices);
    pthread_mutex_destroy(&index_lock);
    SFREE(rdc_argpool);
}

}  /* extern "C */
