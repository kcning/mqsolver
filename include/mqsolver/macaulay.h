/* macaulay.h: header file for structure Macaulay
 */

#ifndef __MACAULAY_H__
#define __MACAULAY_H__

#include <stdint.h>
#include <sys/sysinfo.h>


#include "math.h"
#include "threadpool.h"
#include "rmac.h"

/* ===========================================================================
 *                              structure Macaulay definition
 * ===========================================================================
 */

#define MAC_DROW        (UINT32_MAX)
#define MAC_CROW        (UINT32_MAX-1)

typedef struct {
    union {
        uint64_t* d;
        uint64_t* c;
        uint32_t* s;
    } row;

    // if set to MAC_DROW, then the row is dense. If set to MAC_CROW, then the
    // row is compact sparse. Otherwise the row is sparse and it records the
    // number of terms
    uint32_t term_num;

} MacRow;

typedef struct {
    uint64_t eq_num;
    uint64_t mq_eq_num;
    uint64_t var_num;
    uint64_t term_num;
    uint64_t deg;
    uint32_t* mem_blk;   /* consecutive memory used by sparse rows */
    uint64_t* dmem_blk;  /* consecutive memory used by dense rows */
    uint64_t* cmem_blk;  /* consecutive memory used by compact rows */
    MacRow eqs[];

} Macaulay;

/* ===========================================================================
 *                              function prototypes
 * ===========================================================================
 */

/* function: mac_memsize
 * usage: compute the size of memory block needed for sparse Macaulay
 * arguments:
 *      1) eq_num: number of equations in the original MQ system
 *      2) var_num: number of variables in the original MQ system
 *      3) deg: degree of the Macaulay matrix
 * return: the size of memory needed in bytes
 */
uint64_t
mac_memsize(const uint64_t eq_num, const uint64_t var_num, const uint64_t deg);

/* function: mac_drow_memsize
 * usage: compute how much memory is needed for dense rows
 * argument:
 *      1) term_num: number of terms in a Macaulay
 *      2) drow_num: number of dense rows, which is the number of missing pivots
 * return: memory size in bytes
 */
uint64_t
mac_drow_memsize(const uint64_t term_num, const uint64_t drow_num);

/* function: mac_create
 * usage: create a macaulay matrix container
 * arguments:
 *      1) eq_num: number of equations in the original MQ system
 *      2) var_num: number of variables in the original MQ system
 *      3) deg: degree of the Macaulay matrix
 * return: a pointer to structure Macaulay
 */
Macaulay*
mac_create(const uint64_t eq_num, const uint64_t var_num,
           const uint64_t deg);

/* function: mac_free
 * usage: release a structure Macaulay
 * arguments:
 *      1) mac: pointer to structure Macaulay
 * return: void
 */
void
mac_free(Macaulay* mac);

/* function: mac_row_num
 * usage: given the number of variables and equations in the original MQ
 *      system, compute the number of rows in the Macaulay matrix.
 * arguments:
 *      1) eq_num: number of equations
 *      2) var_num: number of variables
 *      3) deg: degree of the Macaulay matrix
 * return: the number of rows
 */
uint64_t
mac_row_num(const uint64_t eq_num, const uint64_t var_num, const uint64_t deg);

/* function: mac_col_num
 * usage: given the number of variables and degree of the Macaulay matrix,
 *      compute the number of columns in the Macaulay matrix.
 * arguments:
 *      1) var_num: number of variables
 *      2) deg: degree of the Macaulay matrix
 * return: the number of columns
 */
uint64_t
mac_col_num(const uint64_t var_num, const uint64_t deg);

/* function: mac_nl_num
 * usage: given the number of variables to fix, compute the number of
 *      non-linear terms in the remaining variables x1, x2, ... xk that prevent
 *      the remaining system from being linear.
 * arguments:
 *      1) var_num: the number of variables
 *      2) deg: the degree of the Macualay matrix, can be 3 or 4
 *      3) fix_var_num: the number of variable to fix
 * return: the number of non-linear terms in x1 ~ xk
 */
uint64_t
mac_nl_num(const uint64_t var_num, const uint64_t deg, const uint64_t fix_var_num);

/* function: mac_init_perm
 * usage: given a Macaulay matrix container and a MQ system, compute the
 *      Macaulay matrix and permute the columns s.t. non-linear terms are
 *      move to the front and linear terms in the back. The permutation is
 *      stable. The function also records the index of the first non-zero
 *      entry in each row.
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) mqsys: pointer to the input MQ system
 *      3) eq_num: number of equations in the MQ system
 *      4) var_num: number of variables in the MQ system
 *      5) term_num: number of terms in the MQ system, without x_i^2
 *      6) idx_map: an array representing the permutation, idx_map[i] stores
 *              the permuted position for column i
 *      7) clz: container for recording the indices of the first non-zero
 *              entries
 * return: void
 */
void
mac_init_perm(Macaulay* mac, bool* mqsys, const uint64_t eq_num,
              const uint64_t var_num, const uint64_t term_num,
              const uint64_t* const __restrict__ idx_map,
              uint64_t* const __restrict__ clz);

/* function: mac_at
 * usage: return the entry given its position
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) row_idx: row index
 *      3) col_idx: column index
 * return: the bit at the given position
 */
bool
mac_at(Macaulay* mac, const uint64_t row_idx, const uint64_t col_idx);

/* function: mac_zrow
 * usage: check if a row is all zeros
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) row_idx: row index
 * return: true if the row is all zero, otherwise false
 */
bool
mac_zrow(Macaulay* mac, const uint64_t row_idx);

/* function: mac_row_clz
 * usage: count the number of leading zero monomials in a row
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) row_idx: row index
 * return: the number of leading zero monomials
 */
uint64_t
mac_row_clz(Macaulay* mac, const uint64_t row_idx);

/* function: mac_nl_indices
 * usage: given the number of variables to keep, compute the indices of the
 *      non-linear terms
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) indices: container for the indices, must be able to hold at least
 *              the number of non-linear terms
 *      3) kvar_num: number of variables to keep
 * return: void
 */
void
mac_nl_indices(Macaulay* mac, uint64_t* indices, const uint64_t kvar_num);

/* function: mac_sep_indices
 * usage: given the number of variables to keep, separate the indices of the
 *      linear and non-linear terms into two sorted groups
 * arguments:
 *      1) deg: degree of the Macaulay matrix
 *      2) var_num: number of variables
 *      3) nl_indices: container for the indices, must be able to hold at least
 *              the number of non-linear terms
 *      4) l_indices: container for linear terms
 *      5) kvar_num: number of variables to keep
 * return: void
 */
void
mac_sep_indices(const uint64_t deg, const uint64_t var_num,
                uint64_t* const __restrict__ nl_indices,
                uint64_t* const __restrict__ l_indices, const uint64_t kvar_num);

/* function: mac_idx_map
 * usage: given the sorted indices of the linear and non-linear terms, compute
 *      the final position of each column
 * arguments:
 *      1) idx_map: container for the map, must be at least as large as the
 *              number of Macaulay terms
 *      2) nl_indices: indices of non-linear terms
 *      3) l_indices: indices of linear terms
 *      4) nlt_num: number of non-linear terms
 *      5) lt_num: number of linear terms
 * return: void
 */
void
mac_idx_map(uint64_t* const __restrict__ idx_map,
            const uint64_t* const __restrict__ nl_indices,
            const uint64_t* const __restrict__ l_indices,
            const uint64_t nlt_num, const uint64_t lt_num);

/* function: mac_calc_rmac
 * usage: perform row reduction on the left sub-matrix to compute the
 *      sub-system. The Macaulay matrix must be pivoted.
 * arguments:
 *      1) s: pointer to structure RMac
 *      2) mac: pointer to structure Macaulay
 *      3) ltnum: number of linear terms
 *      4) drows: indices of the dense rows
 *      5) drow_num: number of dense pivot rows
 *      6) keq_num: number equations to keep as sub-system candidates
 *      7) tpool: pointer to structure threadpool_t
 * return: void
 */
void
mac_calc_rmac(RMac* s, Macaulay* mac, const uint64_t ltnum,
              const uint32_t* const __restrict__ drows,
              const uint64_t drow_num, const uint64_t keq_num,
              threadpool_t* const tpool);

/* function: mac_pvt
 * usage: perform pivoting on the permuted Macaulay matrix. The routine tries
 *      to move rows to its final position to form a diagonal line before row
 *      reduction. The rows which are not in their final position and the
 *      sub-system candidates are transformed into dense format and their indices
 *      are stored in a container. The last keq_num of of the container are the
 *      indices of sub-system candidates.
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) col_num: number of columns to perform Gaussian elimination
 *      3) keq_num: number equations to reduce below the last pivot
 *      4) clz: the indices of the first non-zero entry of each row
 *      5) drows: container for indices of the dense rows
 *      6) mq_term_ratio: the ratio of terms present in an eq in the original
 *              MQ system
 * return: the number of values in the domain without a suitable row
 */
uint64_t
mac_pvt(Macaulay* __restrict__ mac, const uint64_t col_num,
        const uint64_t keq_num, uint64_t* const __restrict__ clz,
        uint32_t* const __restrict__ drows, const double mq_term_ratio);

/* function: mono_idx
 * usage: given the degree and the variables in a monomial, computes its index
 *      according to grlex order.
 * arguments:
 *      1) deg: degree of the Macaulay matrix
 *      2) var_num: number of variables in the system
 *      3) vars: an sorted array of integers representing the monomial,
 *              e.g. for x1x3x4, one must pass [0, 2, 3]
 *      4) mvar_num: number of variables in the monomial
 * return: the index
 */
uint64_t
mono_idx(const uint64_t deg, const uint64_t var_num,
         const uint64_t* vars, const uint64_t mvar_num);

/* function: fmono_idx
 * usage: given the degree and the multiplier xi (and possibly xj), computes the
 *      index for the first monomials of the product xi (* xj) * eq. If we
 *      organize the rows as block of equation which are multiplied by the same
 *      multipliers, these indices form a boundary in the Macaulay matrix where
 *      all entries below the equation block in the same column and before in the
 *      same row must be zero.
 * arguments:
 *      1) deg: degree of the Macaulay matrix
 *      2) var_num: number of variables
 *      3) vas: an sorted array of integers representing the multipliers, e.g
 *              for x1, x5, one must pass [0, 4]
 *      4) mvar_num: number of multipliers, can be 1 or 2
 * return: the index aforementioned
 */
uint64_t
fmono_idx(const uint64_t deg, const uint64_t var_num,
          const uint64_t* const vars, const uint64_t mvar_num);

/* function: mac_vbound
 * usage: given the index of a row in the Macaulay matrix, return the index of
 *      the first monomial that might be nonzero
 * arguments:
 *      1) eq_num: number of equations in the original MQ system
 *      2) var_num: number of variables in the original MQ system
 *      3) deg: degree of the Macaulay matrix
 *      4) idx: index of the row
 * return: index of the column representing the monomial
 */
uint64_t
mac_vbound(const uint64_t eq_num, const uint64_t var_num, const uint64_t deg,
           const uint64_t idx);

/* function: mac_hsnum
 * usage: compute the number of column where the horizontal boundary changes
 * arguments:
 *      1) var_num: number of variables
 *      2) deg: degree of the Macaulay matrix
 * return: the number
 */
uint64_t
mac_hsnum(const uint64_t var_num, const uint64_t deg);

/* function: mac_hseps
 * usage: compute the column indices where the horizontal boundary changes
 *      and the horizontal boundary
 * arguments:
 *      1) ses: storage for the indices of horizonal separators
 *      2) hbs: storge for the horizonal boundary
 *      3) eq_num: number of eqs in the original MQ system
 *      4) var_num: number of variables in the original MQ system
 *      5) deg: degree of the Macaulay matrix
 * return: void
 */
void
mac_hseps(uint64_t* const __restrict__ ses, uint64_t* const __restrict__ hbs,
          const uint64_t eq_num, const uint64_t var_num, const uint64_t deg);

/* function: mac_hbound
 * usage: given the index of a column in the Macaulay matrix, return the index of
 *      the last equation whose monomial corresponding to that index can be
 *      nonzero.
 * arguments:
 *      1) idx: index of the column
 *      2) senum: size of seps
 *      3) ses: precomputed indices where the horizontal bound changes
 *      4) hbs: horizontal bounds
 * return: the index of the row representing the equation
 */
uint64_t
mac_hbound(const uint64_t idx, const uint64_t senum,
           const uint64_t* const __restrict__ ses,
           const uint64_t* const __restrict__ hbs);

/* function: mac_drow_memsize
 * usage: compute how much memory is needed for dense rows
 * argument:
 *      1) term_num: number of terms in a Macaulay
 *      2) drow_num: number of dense rows, which is the number of missing pivots
 * return: memory size in bytes
 */
uint64_t
mac_drow_memsize(const uint64_t term_num, const uint64_t drow_num);

/* function: mac_crow_memsize
 * usage: compute how much memory is needed for compact rows
 * argument:
 *      1) mq_term_num: number of terms in the original MQ system, without x^2
 *      2) crow_num: number of compact rows, which is the number of pivots
 *              found already
 *      3) mq_term_ratio: the ratio of terms present in an eq in the original
 *              MQ system
 * return: memory size in bytes
 */
uint64_t
mac_crow_memsize(const uint64_t mq_term_num, const uint64_t crow_num,
                 const double mq_term_ratio);

#endif  /* __MACAULAY_H__ */
