/* rmac.h: header file for structure RMac */

#ifndef __RMAC_H__
#define __RMAC_H__


#include <stdint.h>
#include <stdbool.h>

#ifdef MAC_DROW_SSE
    #include <emmintrin.h>
#endif


#include "util.h"
#include "cuda_util.h"
#include "graycode.h"
#include "threadpool.h"

/* the minimal ratio of missing pivot rows for the extracted sub-system to
 * behave like a random system
 */
#define RMAC_MPR_THRESHOLD      (1.0f)

/* ===========================================================================
 *                              structure Macaulay definition
 * ===========================================================================
 */

/* reduced Macaulay matrix */
typedef struct {
    uint64_t eq_num;    /* number of eqs (rows) */
    uint64_t term_num;  /* number of monomials (bits, columns) */
    uint64_t slot_num;  /* number of slots for a row */
    uint64_t* mem;      /* memory block for the rows */
    uint64_t* eqs[];    /* pointer to the start of each row */

} RMac;

/* ===========================================================================
 *                              function prototypes
 * ===========================================================================
 */

/* function: rmac_memsize
 * usage: given the number of eqs and terms, compute the memory needed for rmac
 * arguments:
 *      1) eq_num: number of eqs
 *      2) term_num: number of terms
 * return: size in bytes
 */
__host__ uint64_t
rmac_memsize(const uint64_t eq_num, const uint64_t term_num);

/* function: rmac_create
 * usage: create a container for RMac
 * arguments:
 *      1) eq_num: number of eqs
 *      2) term_num: number of terms
 * return: pointer to the newly created structure RMac
 */
__host__ RMac*
rmac_create(const uint64_t eq_num, const uint64_t term_num);

/* function: rmac_free
 * usage: free structure RMac
 * arguments:
 *      1) sys: pointer to structure RMac
 * return: void
 */
__host__ void
rmac_free(RMac* sys);

/* function: rmac_row
 * usage: given the row index, return the starting address for the row in RMac
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) idx: row index
 * return: the address
 */
__host__ uint64_t*
rmac_row(RMac* const sys, const uint64_t idx);

/* function: rmac_at
 * usage: given the row and column index, return the entry
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) ridx: row index
 *      3) cidx: column index
 * return: the entry
 */
__host__ bool
rmac_at(RMac* const sys, const uint64_t ridx, const uint64_t cidx);

/* function: rmac_gauss_elim
 * usage: perform Gaussian elimination on the sub-system
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) col_num: number of columns to perform elimination
 * return: void
 */
__host__ void
rmac_gauss_elim(RMac* const sys, const uint64_t col_num);

/* function: rmac_zrow
 * usage: check if a row is all zeros
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) row_idx: row index
 * return: true if the row is all zero, otherwise false
 */
__host__ bool
rmac_zrow(RMac* sys, const uint64_t row_idx);

/* function: rmac_perm_rows
 * usage: permute the rows according to an index map
 * arguments:
 *      1) sys: pointer to structure RMac
 *      2) idx_map: row index map
 * return: void
 */
__host__ void
rmac_perm_rows(RMac* sys, const uint32_t* const idx_map);

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
              uint32_t* const __restrict__ host_rcount);

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
              threadpool_t* const tpool);

#endif  /* __RMAC_H__ */
