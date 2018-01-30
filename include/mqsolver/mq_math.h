/* mq_math.h: define math functions */

#ifndef __MQ_MATH_H__
#define __MQ_MATH_H__

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>  /* for debug */


#include "util.h"

/* binomial coefficients */
#define binom2(n) \
    ( ((n) * ((n)-1)) / 2)

#define binom3(n) \
    ( ( ( ((n) * ((n)-1)) / 2 ) * ((n)-2) ) / 3 )

/* Note that integer overflow might happen when n > 400. For the target
 * problem instances, n << 400 so it should be fine.
 */
#define binom4(n) \
    ( ( ( ( ( ((n) * ((n)-1)) / 2 ) * ((n)-2) ) / 3 ) * ((n)-3) ) / 4)

/* macros: midx[0-3]
 * usage: given the degree and the variables in a monomial, computes its index
 *      according to grlex order.
 * arguments:
 *      1) deg: degree of the Macaulay matrix
 *      2) var_num: number of variables in the system
 *      3-5) vars: sorted integers representing the monomial,
 *              e.g. for x1x3x4, one must pass 0, 2, 3
 * return: the index
 */
#define midx0(deg, vnum) \
    (((3 == (deg)) ? 0 : binom4(vnum)) + binom3(vnum) + binom2(vnum) + (vnum))

#define midx1(deg, vnum, var_idx) \
    (((3 == (deg)) ? 0 : binom4(vnum)) + binom3(vnum) + binom2(vnum) + (var_idx))

#define midx2(deg, vnum, var1_idx, var2_idx) \
    (((3 == (deg)) ? 0 : binom4(vnum)) + binom3(vnum) + binom2(var2_idx) + (var1_idx))

#define midx3(deg, vnum, var1_idx, var2_idx, var3_idx) \
    (((3 == (deg)) ? 0 : binom4(vnum)) + binom3(var3_idx) + binom2(var2_idx) + (var1_idx))

/* macros: midx4
 * usage: given the variables in a degree 4 monomial, computes its index
 *      according to grlex order.
 * arguments:
 *      1) var_num: number of variables in the system
 *      2-5) vars: sorted integers representing the monomial,
 *              e.g. for x1x3x4x6, one must pass 0, 2, 3, 5
 * return: the index
 */
#define midx4(vnum, var1_idx, var2_idx, var3_idx, var4_idx) \
    (binom4(var4_idx) + binom3(var3_idx) + binom2(var2_idx) + (var1_idx))

/* function: binom
 * usage: compute binominal coefficient
 * arguments:
 *      1) n: the size of the set
 *      2) k: the size of the subset
 * return: the binominal coefficient C(n, k)
 */
uint64_t
binom(uint32_t n, uint32_t k);

/* function: sum_binom
 * usage: compute the sum of binominal coefficients
 * arguments:
 *      1) n: the size of the set
 *      2) k: the size of the subset
 * return: the sum of binominal coefficients C(n, i), for i in 0...k
 */
uint64_t
sum_binom(uint32_t n, uint32_t k);

/* function: diff_eq
 * usage: differentiate the given polynomial with respect to one varible
 *      in GF(2)
 * arguments:
 *      1) func: pointer to a bool array, which are the coefficients
 *              of the polynomial in graded reverse lex order
 *      2) term_num: number of the terms
 *      3) var_num: number of variables
 *      4) idx: which variable x_{i+1} to differentiate against. 
 *      5) result: a pointer to a bool array which will hold the result.
 *              The size must be at least var_num + 1
 * return: void
 */
void
diff_eq(bool* func, uint64_t term_num, uint64_t var_num, uint64_t idx,
        bool* result);

/* function: sum_inc
 * usage: compute \sigma_0^n
 * argument: upperbound of the summation
 * return: the sum
 */

uint64_t
sum_inc(uint64_t n);

/* function: sum_dec
 * usage: summation of decreasing non-negative integers, e.g. 3, 2, 1, 0
 * arguments:
 *      1) v: initial value
 *      2) n: number of values to sum in the sequence
 * return: sum of the sequence
 */
uint64_t
sum_dec(uint64_t v, uint64_t n);

/* function: gj_elim_trans
 * usage: given a transposed linear system, reduce it to reduced row echelon
 *      form using Gauss-Jordan elimination. The system is stored as an array
 *      of uint64_t where each uint64_t holds the same column for all rows.
 * arguments:
 *      1) sys: the transposed linear system
 *      2) col_num: number of columns
 * return: void
 */
void
gj_elim_trans(uint64_t* sys, const uint64_t col_num);

/* function: solve_ls_trans
 * usage: given a transposed linear system, perform Gaussian-Jordan elimination
 *      to solve the linear system. The system is stored as an array of
 *      uint64_t where each uint64_t holds the column for all rows. Note that
 *      the system is modified in place.
 * arguments:
 *      1) sys: the transposed linear system
 *      2) col_num: number of columns
 * return: true if the system has a unique or multiple solutions, false
 *      otherwise
 */
bool
solve_ls_trans(uint64_t* sys, const uint64_t col_num);

/* function: check_ls_trans
 * usage: given a transposed linear system, check if the system is solvable
 *      The system is stored as an array of uint64_t where each uint64_t holds
 *      the column for all rows. The function will abort early if the system is
 *      found to be inconsistent.
 * arguments:
 *      1) sol: a temporary container, should be at least
 *              col_num * uint64_t * 8 bytes.
 *      2) sys: the transposed linear system
 *      3) col_num: number of columns
 *      4) eq_num: number of equations
 * return: true if the system has a unique or multiple solutions, false
 *      otherwise
 */
bool
check_ls_trans(uint64_t* sol, const uint64_t* const sys, const uint64_t col_num,
               const uint64_t eq_num);

/* function: det
 * usage: given a transposed lienar system in its reduced row echelon form,
 *      compute its determinant
 * arguments:
 *      1) sys: the transposed linear system
 *      2) col_num: number of columns
 * return: the determinant
 */
bool
det(uint64_t* sys, const uint64_t col_num);

/* function: dpcount
 * usage: given a matrix stored in column-wise foramt, perform pop count on
 *      the diagonal line of a matrix
 * arguments:
 *      1) sys: the matrix stored in column-wise format, where each uint64_t
 *              holds the same column for all the rows
 *      2) col_num: number of columns
 * return: the number of diagonal entries that are set to 1
 */
uint64_t
dpcount(const uint64_t* const sys, const uint64_t col_num);

/* ===========================================================================
 *                              generic macros definition
 * ===========================================================================
 */

/* macro: max
 * usage: find the maximum value in a dataset
 *      1) stats: an array of values. can be of type uint8_t, double
 *      2) num: the size of the array
 * retrun: the maximum, same type as elements in stats
 */
#define max(stats, num) _Generic((stats), \
    uint8_t* : max_i8, \
    double*  : max_f64 \
    )(stats, num)

/* macro: avg
 * usage: compute the average of a dataset
 * arguments:
 *      1) stats: an array of values. can be of type uint8_t, double
 *      2) num: the size of the array
 * return: the average as a double value
 */
#define avg(stats, num) _Generic((stats), \
    uint8_t* : avg_i8, \
    double*  : avg_f64 \
    )(stats, num)

/* macro: std_dev
 * usage: compute the standard deviation of a dataset
 * arguments:
 *      1) stats: an array of values, can be of type uint8_t, double
 *      2) num: the size of the array
 *      3) avg: the average value of the dataset
 * return: the stddev as a double value
 */
#define std_dev(stats, num, average) _Generic((stats), \
    uint8_t* : std_dev_i8, \
    double*  : std_dev_f64 \
    )(stats, num, average)

/* ===========================================================================
 *                              boilerplates for generic macros
 *                              no need to call those directly
 * ===========================================================================
 */

/* function: sum_i8
 * usage: compute the sum of a dataset
 * arguments:
 *      1) stats: an array of uint8_t values
 *      2) num: the size of the array
 * return: the sum
 */
uint64_t
sum_i8(uint8_t* stats, const uint64_t num);

/* function: avg_f64
 * usage: compute the average of a dataset
 * arguments:
 *      1) stats: an array of double values
 *      2) num: the size of the array
 * return: the average
 */
double
avg_f64(double* stats, const uint64_t num);

/* function: avg_i8
 * usage: compute the average of a dataset
 * arguments:
 *      1) stats: an array of uint8_t values
 *      2) num: the size of the array
 * return: the average
 */
double
avg_i8(uint8_t* stats, const uint64_t num);

/* function: std_dev_f64
 * usage: compute the standard deviation of a dataset
 * arguments:
 *      1) stats: an array of double values
 *      2) num: the size of the array
 *      3) avg: the average value of the dataset
 * return: the stddev
 */
double
std_dev_f64(double* stats, const uint64_t num, double average);

/* function: std_dev_i8
 * usage: compute the standard deviation of a dataset
 * arguments:
 *      1) stats: an array of uint8_t values
 *      2) num: the size of the array
 *      3) avg: the average value of the dataset
 * return: the stddev
 */
double
std_dev_i8(uint8_t* stats, const uint64_t num, double average);

/* function: max_i8
 * usage: find the maximum value in a dataset
 *      1) stats: an array of uint8_t values
 *      2) num: the size of the array
 * retrun: the maximum
 */
uint8_t
max_i8(uint8_t* stats, const uint64_t num);

/* function: max_f64
 * usage: find the maximum value in a dataset
 *      1) stats: an array of double values
 *      2) num: the size of the array
 * retrun: the maximum
 */
double
max_f64(double* stats, const uint64_t num);

#endif  /* __MQ_MATH_H__ */
