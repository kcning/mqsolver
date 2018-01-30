/* debug.c: implementation of debug.h */

#include "debug.h"

/* function: gaussian_elim
 * usage: given a linear system, reduce it to row echelon form using Gaussian
 *      elimination. Note the system in modified in place.
 *      1) sys: the linear system
 *      2) row_num: number of rows
 *      3) col_num: number of columns
 * return: void
 */
inline void
gaussian_elim(bool* sys, uint64_t row_num, uint64_t col_num) {
    reduce_lsub(sys, row_num, col_num, col_num);
}

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
reduce_lsub(bool* sys, uint64_t row_num, uint64_t col_num, uint64_t lcol_num) {
    assert(lcol_num <= col_num);
    const uint64_t bound = (row_num < lcol_num) ? row_num : lcol_num;
    uint64_t i, j, k;
    bool tmp_row[col_num];
    for(i = 0; i < bound; ++i) {
        // find the i-th pivot
        for(j = i; j < row_num; ++j) {
            if(true == sys[j * col_num + i]) { // sys[j][i]
                break;
            }
        }
        if(row_num == j) { // singular
            continue;
        }

        // swap i-th and j-th rows
        memcpy(tmp_row, sys + i * col_num, col_num);
        memcpy(sys + i * col_num, sys + j * col_num, col_num);
        memcpy(sys + j * col_num, tmp_row, col_num);

        // for all the rows below pivot
        for(j = i+1; j < row_num; ++j) {
            if(true == sys[j * col_num + i]) { // sys[j][i]
                // subtract i-th row from the row
                for(k = 0; k < col_num; ++k) {
                    sys[j * col_num + k] ^= sys[i * col_num + k];
                }
            }
        }
    }
}

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
check_solvable(bool* sys, uint64_t row_num, uint64_t col_num) {
    assert(row_num >= col_num);
    gauss_jordan_elim(sys, row_num, col_num);
    return !sys[ (col_num-1) * col_num + col_num-1];
}

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
rref_lsub(bool* sys, uint64_t row_num, uint64_t col_num, uint64_t lcol_num) {
    assert(lcol_num <= col_num);
    const uint64_t bound = (row_num < lcol_num) ? row_num : lcol_num;
    uint64_t i, j, k;
    bool tmp_row[col_num];
    for(i = 0; i < bound; ++i) {
        // find the i-th pivot
        for(j = i; j < row_num; ++j) {
            if(true == sys[j * col_num + i]) { // sys[j][i]
                break;
            }
        }
        if(row_num == j) { // singular
            continue;
        }

        // swap i-th and j-th rows
        memcpy(tmp_row, sys + i * col_num, col_num);
        memcpy(sys + i * col_num, sys + j * col_num, col_num);
        memcpy(sys + j * col_num, tmp_row, col_num);

        // for all the rows below pivot
        for(j = i+1; j < row_num; ++j) {
            if(true == sys[j * col_num + i]) { // sys[j][i]
                // subtract i-th row from the row
                for(k = 0; k < col_num; ++k) {
                    sys[j * col_num + k] ^= sys[i * col_num + k];
                }
            }
        }

        // for all the rows above pivot
        for(j = 0; j < i; ++j) {
            if(true == sys[j * col_num + i]) { // sys[j][i]
                // subtract i-th row from the row
                for(k = 0; k < col_num; ++k) {
                    sys[j * col_num + k] ^= sys[i * col_num + k];
                }
            }
        }
    }
}

/* function: gauss_jordan_elim
 * usage: given a linear system, reduce it to reduced row echelon form using
 *      Gauss-Jordan elimination. Note the system in modified in place.
 * arguments:
 *      1) sys: the linear system
 *      2) row_num: number of rows
 *      3) col_num: number of columns
 * return: void
 */
inline void
gauss_jordan_elim(bool* sys, uint64_t row_num, uint64_t col_num) {
    rref_lsub(sys, row_num, col_num, col_num);
}

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
              const uint64_t deg) {
    for(uint64_t i = 0; i < lsys_eq_num; ++i) {
        // constant and linear terms
        memcpy(lsys + i * (var_num+1),
               mac + (i+1) * mac_term_num - (var_num+1),
               var_num+1);

        // fix quadratic terms
        uint64_t deg2_idx = (i+1) * mac_term_num - (var_num+1);
        for(uint64_t j = 0; j < fix_var_num; ++j) {
            deg2_idx -= var_num-j-1;
            if(count & (0x1UL << j)) { // the var is set to 1
                // xor the resultant linear terms by setting var to 1
                for(uint64_t k = 0; k < var_num-1-j; ++k) {
                    lsys[i * (var_num+1) + k] ^= mac[deg2_idx + k];
                }
            }
        }

        // fix degree 3 monomials
        uint64_t deg3_idx = (i+1) * mac_term_num - (var_num+1) - binom(var_num, 2);
        for(uint64_t j = 0; j < fix_var_num; ++j) { // last var
            if(count & 0x1UL << j) { // the var is set to 1
                for(uint64_t k = 0; k < fix_var_num-1-j; ++k) { // 2nd last var
                    deg3_idx -= var_num-2-j-k;
                    if(count & 0x1UL << (k+1+j)) { // the 2nd var is set to 1
                        for(uint64_t l = 0; l < var_num-2-j-k; ++l) { // first var
                            // xor linear term
                            lsys[i * (var_num+1) + l] ^= mac[deg3_idx + l];
                                
                        }
                    }
                }

                for(uint64_t k = fix_var_num-1-j; k < var_num-1-j; ++k) {
                    deg3_idx -= var_num-2-j-k;
                }
            } else {
                deg3_idx -= binom(var_num-1-j, 2);
            }
        }

        // fix degree 4 monomials
        if(4 == deg) {
            uint64_t deg4_idx = (i+1) * mac_term_num - 1 - var_num
                - binom(var_num, 2) - binom(var_num, 3) - 1;
            for(uint64_t j = 0; j < fix_var_num; ++j) { // last var in the monomial
                // 2nd last var in the monomial
                for(uint64_t k = 0; k < fix_var_num-1-j; ++k) {
                    // 3rd last var in the monomial
                    for(uint64_t l = 0; l < fix_var_num-2-j-k; ++l) {
                        // can be reduced to linear term
                        for(uint64_t m = 0; m < var_num-3-j-k-l; ++m) {
                            lsys[i * (var_num+1) + var_num-3-j-k-l-m-1] ^= 
                                (count >> j) & (count >> (1+j+k)) &
                                (count >> (2+j+k+l)) & mac[deg4_idx-m];
                        }
                        // move deg4_idx to the last monomial of the form
                        // x_?x_{var_num-2-j-k-l-1}x_{var_num-1-j-k}x_{var_num-j}
                        deg4_idx -= var_num-3-j-k-l;
                    }

                    for(uint64_t l = fix_var_num-2-j-k; l < var_num-3-j-k; ++l) {
                        deg4_idx -= var_num-3-j-k-l;
                    }
                }

                for(uint64_t k = fix_var_num-1-j; k < var_num-3-j; ++k) {
                    for(uint64_t l = 0; l < var_num-3-j-k; ++l) {
                        deg4_idx -=  var_num-3-j-k-l;
                    }
                }
            }
        }


        // finally fix linear terms
        for(uint64_t j = 0; j < fix_var_num; ++j) {
            if(true == lsys[i * (var_num+1) + (var_num-1-j)] &&
              (count & (0x1UL << j))) {
                // flip constant term
                lsys[i * (var_num+1) + var_num] ^= true;
            }
            // set the fixed linear term to zero
            lsys[i * (var_num+1) + (var_num-1-j)] = false;
        }
    }
}
