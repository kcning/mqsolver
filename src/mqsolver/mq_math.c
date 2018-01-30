/* mq_math.c: implementation of mq_math.h */

#include "mq_math.h"

/* function: binom
 * usage: compute binominal coefficient
 * arguments:
 *      1) n: the size of the set
 *      2) k: the size of the subset
 * return: the binominal coefficient C(n, k)
 */
uint64_t
binom(uint32_t n, uint32_t k) {
    if(n < k) {
        return 0;
    }

    if(n == k || 0 == k) {
        return 1;
    }

    if((n-k) < k) k = n - k;

    uint64_t res = 1;
    uint32_t i = 0;
    for(i = 0; i < k; ++i) {
        res *= n-i;
        res /= i+1;
    }

    return res;
}

/* function: sum_binom
 * usage: compute the sum of binominal coefficients
 * arguments:
 *      1) n: the size of the set
 *      2) k: the size of the subset
 * return: the sum of binominal coefficients C(n, i), for i in 0...k
 */
inline uint64_t
sum_binom(uint32_t n, uint32_t k) {
    uint64_t sum = 0;

    uint32_t i = 0;
    for(i = 0; i <= k; ++i) {
        sum += binom(n, i);
    }

    return sum;
}

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
        bool* result) {
    assert(idx < var_num);

    // set all terms to zero
    memset(result, 0x0, sizeof(bool) * (var_num + 1));

    // ignore constant term in func, it's gone

    // diff x_{i+1}, which becomes the constant term of the result
    result[var_num] = func[term_num - 1 - (var_num - idx)];

    // diff x1x2, x2x3, ...
    uint64_t cursor = term_num - 1 - var_num - 1;  // start with x_{var_num}^2
    uint64_t i, j;
    uint64_t bound = (0 == idx) ? 1 : idx;
    for(i = var_num - 1; i >= bound; --i) {  // check x_?x_{i+1}, and skip x1^2
        if(i == idx) {
            for(j = 1; j <= i; ++j) {
                result[i - j] ^= func[cursor - j];
            }
        } else {
            result[i] ^= func[cursor - (i - idx)];
        }
        cursor -= i+1;  // move to x_{i-1}^2
    }
}

/* function: sum_inc
 * usage: compute \sigma_0^n
 * argument: upperbound of the summation
 * return: the sum
 */

inline uint64_t
sum_inc(uint64_t n) {
    uint64_t i;
    uint64_t sum = 0;
    for(i = 1; i <= n; ++i) {
        sum += i;
    }
    return sum;
}

/* function: sum_dec
 * usage: summation of decreasing non-negative integers, e.g. 3, 2, 1, 0
 * arguments:
 *      1) v: initial value
 *      2) n: number of values to sum in the sequence
 * return: sum of the sequence
 */
uint64_t
sum_dec(uint64_t v, uint64_t n) {
    uint64_t i;
    uint64_t sum = 0;
    for(i = 0; i < n; ++i) {
        sum += v;
        --v;
    }
    return sum;
}

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
gj_elim_trans(uint64_t* sys, const uint64_t col_num) {
    for(uint64_t i = 0; i < col_num; ++i) {
        
        // NOTE: one can use this flag to avoid branching
        //bool sf = sys[i] >> i; // skip flag
        //uint64_t skip_mask = (sf ^ 0x1UL) - 0x1UL; // sign extend sf
        // all 1's if don't skip, 0x0UL if skip

        if( !(sys[i] >> i) ) {
            continue; // singular
        }

        // find pivot position, which is the number of trailing zeros
        const int piv = __builtin_ctzll(sys[i] >> i) + i;

        //uint64_t mask = (sys[i] ^ (0x1UL << piv)) & skip_mask;
        const uint64_t mask = sys[i] ^ (0x1UL << piv);
        sys[i] ^= mask;

        // row reduction
        for(uint64_t j = i+1; j < col_num; ++j) {
            uint64_t b = (sys[j] >> piv) & 0x1UL;
            //sys[j] ^= mask & (~b + 1);
            sys[j] ^= mask & ((b ^ 0x1UL) - 0x1UL); // sign extend b
        }

        // swap rows, i.e. swap all the ith bits with piv-th bits
        for(uint64_t j = 0; j < col_num; ++j) {
            // i-th bit ^ piv-th bit
            uint64_t tmp = ((sys[j] >> i) ^ (sys[j] >> piv)) & 0x1UL;
            //sys[j] ^= ((tmp << i) | (tmp << piv)) & skip_mask;
            sys[j] ^= (tmp << i) | (tmp << piv);
        }
    }
}

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
inline bool
solve_ls_trans(uint64_t* sys, const uint64_t col_num) {
    gj_elim_trans(sys, col_num);
    return !((sys[col_num-1] >> (col_num-1)) & 0x1UL);
}

/* subroutine of check_ls_trans, extract an equation from the transposed
 * linear system
 */
static inline uint64_t
extract_eq(const uint64_t* const sys, const uint64_t eq_idx, const uint64_t col_num) {
    uint64_t eq = 0x0UL;

    // format of eq: 0000...01x_kx_{k-1}x_...x_2x_1
    for(uint64_t i = 0; i < col_num; ++i) {
        eq |= ((sys[i] >> eq_idx) & 0x1UL) << i;
    }

    return eq;
}

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
               const uint64_t eq_num) {

    memset(sol, 0x0, sizeof(uint64_t) * col_num); // reset
    const uint64_t ctr = 0x1UL << (col_num-1);
    uint64_t ieq_num = 0; // number of linear independent eqs

    for(uint64_t eq_idx = 0; eq_idx < eq_num; ++eq_idx) {
        uint64_t eq = extract_eq(sys, eq_idx, col_num);
        while(eq) {
            if(eq == ctr) { // 0 = 1
                return false;
            }

            int pos = __builtin_ctzll(eq);
            if(sol[pos]) {
                eq ^= sol[pos];
            } else {
                sol[pos] = eq;
                ++ieq_num;

                eq = 0x0UL;
            }
        }
    }

    // multiple solutions
    return true;
}

/* function: det
 * usage: given a transposed lienar system in its reduced row echelon form,
 *      compute its determinant
 * arguments:
 *      1) sys: the transposed linear system
 *      2) col_num: number of columns
 * return: the determinant
 */
inline bool
det(uint64_t* sys, const uint64_t col_num) {
    uint64_t d = 0x1UL;

    for(uint64_t i = 0; i < col_num; ++i) {
        d &= sys[i] >> i;
    }

    return d;
}

/* function: dpcount
 * usage: given a matrix stored in column-wise foramt, perform pop count on
 *      the diagonal line of a matrix
 * arguments:
 *      1) sys: the matrix stored in column-wise format, where each uint64_t
 *              holds the same column for all the rows
 *      2) col_num: number of columns
 * return: the number of diagonal entries that are set to 1
 */
inline uint64_t
dpcount(const uint64_t* const sys, const uint64_t col_num) {
    uint64_t c = 0;

    for(uint64_t i = 0; i < col_num; ++i) {
        c += (sys[i] >> i) & 0x1UL;
    }

    return c;
}

/* function: sum_i8
 * usage: compute the sum of a dataset
 * arguments:
 *      1) stats: an array of uint8_t values
 *      2) num: the size of the array
 * return: the sum
 */
inline uint64_t
sum_i8(uint8_t* stats, const uint64_t num) {
    uint64_t s = 0;

    for(uint64_t i = 0; i < num; ++i) {
        s += stats[i];
    }

    return s;
}

/* function: avg_f64
 * usage: compute the average of a dataset
 * arguments:
 *      1) stats: an array of double values
 *      2) num: the size of the array
 * return: the average
 */
inline double
avg_f64(double* stats, const uint64_t num) {
    double average = 0.0;

    for(uint64_t i = 0; i < num; ++i) {
        average += stats[i];
    }

    return average / num;
}

/* function: avg_i8
 * usage: compute the average of a dataset
 * arguments:
 *      1) stats: an array of uint8_t values
 *      2) num: the size of the array
 * return: the average
 */
inline double
avg_i8(uint8_t* stats, const uint64_t num) {
    return ((double) sum_i8(stats, num)) / num;
}

/* function: std_dev_f64
 * usage: compute the standard deviation of a dataset
 * arguments:
 *      1) stats: an array of double values
 *      2) num: the size of the array
 *      3) avg: the average value of the dataset
 * return: the stddev
 */
inline double
std_dev_f64(double* stats, const uint64_t num, double average) {
    double sdev = 0.0;

    for(uint64_t i = 0; i < num; ++i) {
        sdev += pow(stats[i] - average, 2);
    }

    return sqrt(sdev);
}

/* function: std_dev_i8
 * usage: compute the standard deviation of a dataset
 * arguments:
 *      1) stats: an array of uint8_t values
 *      2) num: the size of the array
 *      3) avg: the average value of the dataset
 * return: the stddev
 */
inline double
std_dev_i8(uint8_t* stats, const uint64_t num, double average) {
    double sdev = 0.0;

    for(uint64_t i = 0; i < num; ++i) {
        sdev += pow( (double) stats[i] - average, 2);
    }

    return sqrt(sdev);
}

/* function: max_i8
 * usage: find the maximum value in a dataset
 *      1) stats: an array of uint8_t values
 *      2) num: the size of the array
 * retrun: the maximum
 */
inline uint8_t
max_i8(uint8_t* stats, const uint64_t num) {
    uint8_t max = 0;

    for(uint64_t i = 0; i < num; ++i) {
        if(max < stats[i]) {
            max = stats[i];
        }
    }

    return max;
}

/* function: max_f64
 * usage: find the maximum value in a dataset
 *      1) stats: an array of double values
 *      2) num: the size of the array
 * retrun: the maximum
 */
inline double
max_f64(double* stats, const uint64_t num) {
    double max = 0;

    for(uint64_t i = 0; i < num; ++i) {
        if(max < stats[i]) {
            max = stats[i];
        }
    }

    return max;
}
