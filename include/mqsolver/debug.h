/* debug.h: some naive implementation of math functions for debugging */

#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <stdbool.h>
#include <stdint.h>
#include <assert.h>  /* for debug */

#include "mq_math.h"
#include "util.h"

/* function: gaussian_elim
 * usage: given a linear system, reduce it to row echelon form using Gaussian
 *      elimination. Note the system in modified in place.
 *      1) sys: the linear system
 *      2) row_num: number of rows
 *      3) col_num: number of columns
 * return: void
 */
void
gaussian_elim(bool* sys, uint64_t row_num, uint64_t col_num);

/* function: reduce_lsub
 * usage: given a linear system, reduce its left submatrix into its row echelon
 *      form using Gaussian elimination. Note the system is modified in place.
 *      1) sys: the linear system
 *      2) row_num: number of rows
 *      3) col_num: number of columns
 *      4) lcol_num: number of columns in the left sub matrix
 * return: void
 */
void
reduce_lsub(bool* sys, uint64_t row_num, uint64_t col_num, uint64_t lcol_num);

/* function: check_solvable
 * usage: given a linear system, check if the system is solvable. Note the
 *      system is modified in place.
 * arguments:
 *      1) sys: the linear system
 *      2) row_num: number of rows
 *      3) col_num: number of columns
 * return: true if solvable, otherwise false
 */
bool
check_solvable(bool* sys, uint64_t row_num, uint64_t col_num);

/* function: rref_lsub
 * usage: given a linear system, try to reduce its left submatrix into its
 *      reduced row echelon form using Gauss-Jordan elimination. Note the
 *      system is modified in place.
 * arguments:
 *      1) sys: the linear system
 *      2) row_num: number of rows
 *      3) col_num: number of columns
 *      4) lcol_num: number of columns in the left sub matrix
 * return: void
 */
void
rref_lsub(bool* sys, uint64_t row_num, uint64_t col_num, uint64_t lcol_num);

/* function: gauss_jordan_elim
 * usage: given a linear system, reduce it to reduced row echelon form using
 *      Gauss-Jordan elimination. Note the system in modified in place.
 * arguments:
 *      1) sys: the linear system
 *      2) row_num: number of rows
 *      3) col_num: number of columns
 * return: void
 */
void
gauss_jordan_elim(bool* sys, uint64_t row_num, uint64_t col_num);

/* function: fix_last_vars
 * usage: fix the last few variables in an equation and turn the input system
 *      into a linear system
 * arguments:
 *      1) lsys: storage for the result
 *      2) mac: input system
 *      3) lsys_eq_num: number of equation for the input system
 *      4) mac_term_num: number of terms in the input system
 *      5) var_num: number of variables
 *      6) fix_var_num: number of variables to fix
 *      7) count: an integer where the last bit represents the last variable
 *              to fix(x_n). For the rest of the bits, idem.
 *      8) deg: degree of the input system
 * return: void
 */
void
fix_last_vars(bool* lsys, bool* mac, const uint64_t lsys_eq_num,
              const uint64_t mac_term_num, const uint64_t var_num,
              const uint64_t fix_var_num, const uint64_t count,
              const uint64_t deg);

#endif  /* __DEBUG_H__ */
