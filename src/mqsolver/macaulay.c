/* macaulay.c: implementation of structure Macaulay
 */

#include "macaulay.h"
#include "drow.h"

/* function: mac_memsize
 * usage: compute the size of memory block needed for sparse Macaulay
 * arguments:
 *      1) eq_num: number of equations in the original MQ system
 *      2) var_num: number of variables in the original MQ system
 *      3) deg: degree of the Macaulay matrix
 * return: the size of memory needed in bytes
 */
inline uint64_t
mac_memsize(const uint64_t eq_num, const uint64_t var_num, const uint64_t deg) {
    assert(deg == 3 || 4 == deg);

    const uint64_t mac_eq_num = mac_row_num(eq_num, var_num, deg);
    const uint64_t mq_term_num = binom(var_num, 2) + var_num + 1;
    uint64_t memsize = sizeof(Macaulay) + sizeof(MacRow) * mac_eq_num;
    return memsize + mac_eq_num * mq_term_num * sizeof(uint32_t);
}

/* function: mac_create
 * usage: create a macaulay matrix container
 * arguments:
 *      1) eq_num: number of equations in the original MQ system
 *      2) var_num: number of variables in the original MQ system
 *      3) deg: degree of the Macaulay matrix
 * return: a pointer to structure Macaulay
 */
Macaulay*
mac_create(const uint64_t eq_num, const uint64_t var_num, const uint64_t deg) {
    assert(3 == deg || 4 == deg);

    const uint64_t mac_eq_num = mac_row_num(eq_num, var_num, deg);
    Macaulay* mac = (Macaulay*) malloc(sizeof(Macaulay) +
                                       sizeof(MacRow) * mac_eq_num);
    if(NULL == mac) {
        return NULL;
    }

    mac->mq_eq_num = eq_num;
    mac->eq_num = mac_eq_num;
    mac->var_num = var_num;
    mac->deg = deg;
    mac->term_num = mac_col_num(var_num, deg);

    const uint64_t mq_term_num = sum_binom(var_num, 2);
    mac->mem_blk = (uint32_t*) malloc(mac_eq_num * mq_term_num * sizeof(uint32_t));
    mac->dmem_blk = NULL;
    mac->cmem_blk = NULL;

    for(uint64_t i = 0; i < mac_eq_num; ++i) {
        mac->eqs[i].term_num = 0;
        mac->eqs[i].row.s = mac->mem_blk + i * mq_term_num;
    }

    return mac;
}

/* function: mac_free
 * usage: release a structure Macaulay
 * arguments:
 *      1) mac: pointer to structure Macaulay
 * return: void
 */
inline void
mac_free(Macaulay* mac) {
    if(mac->dmem_blk) {
        free(mac->dmem_blk);
    }

    if(mac->cmem_blk) {
        free(mac->cmem_blk);
    }

    if(mac->mem_blk) {
        free(mac->mem_blk);
    }
    free(mac);
}

/* function: mac_drow_memsize
 * usage: compute how much memory is needed for dense rows
 * argument:
 *      1) term_num: number of terms in a Macaulay
 *      2) drow_num: number of dense rows, which is the number of missing pivots
 * return: memory size in bytes
 */
inline uint64_t
mac_drow_memsize(const uint64_t term_num, const uint64_t drow_num) {
    return drow_num * sizeof(uint64_t) * drow_slots_num(term_num);
}

/* function: mac_to_drow
 * usage: transform a row from sparse format into dense format
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) idx: index of the row
 *      3) mem: memory block for the dense row, should be all zero
 * return: void
 */
static inline void
mac_to_drow(Macaulay* mac, const uint64_t idx, uint64_t* const mem) {

    for(uint64_t i = 0; i < mac->eqs[idx].term_num; ++i) {
        drow_set_at(mem, mac->eqs[idx].row.s[i]);
    }

    mac->eqs[idx].term_num = MAC_DROW;
    mac->eqs[idx].row.d = mem;
}

/* ===========================================================================
 *                              Compact row
 * ===========================================================================
 */

/* function: crow_slot_num
 * usage: given the number of terms in MQ system, compute the number of
 *      uint64_t needed for a compact row
 * arguments:
 *      1) mq_tnum: the number of terms in original MQ system, without x^2
 *      2) mq_term_ratio: the ratio of terms present in an eq in the original
 *              MQ system
 * return: slot number
 */
static inline uint64_t
crow_slot_num(const uint64_t mq_tnum, const double mq_term_ratio) {
    return 1 + (uint64_t)(mq_tnum * 2 * mq_term_ratio);
}

/* function: crow_idx
 * usage: given the slot idx, return the address to the index for the
 *      monomials in a drow
 * arguments:
 *      1) crow: pointer to the compact row
 *      2) idx: slot index
 * return: address to the slot index where the 64 monomials should be
 *      stored in a drow
 */
static inline uint64_t*
crow_idx(uint64_t* const crow, const uint64_t idx) {
    return crow + 1 + 2 * idx;
}

/* function: crow_mono
 * usage: given the slot idx, address to the monomials stored in the slot
 * arguments:
 *      1) crow: pointer to the compact row
 *      2) idx: slot index
 * return: address to the uint64_t representing the 64 monomials stored
 *      in the slot
 */
static inline uint64_t*
crow_mono(uint64_t* const crow, const uint64_t idx) {
    return crow + 1 + 2 * idx + 1;
}

/* function: crow_snum
 * usage: return the address to the number of slots in a compact row
 * arguments:
 *      1) crow: pointer to the compact row
 * return: address to the number of slots in the crow
 */
static inline uint64_t*
crow_snum(uint64_t* const crow) {
    return crow;
}

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
inline uint64_t
mac_crow_memsize(const uint64_t mq_term_num, const uint64_t crow_num,
                 const double mq_term_ratio) {
    return crow_num * crow_slot_num(mq_term_num, mq_term_ratio) * sizeof(uint64_t);
}

/* function: mac_to_crow
 * usage: transform a row from sparse format into compact format
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) idx: index of the row
 *      3) mem: memory block for the compact row
 *      4) cr_slot_num: max number of slots in the crow
 *      4) drow: a drow
 *      5) dr_slot_num: number of slots in a drow
 * return: void
 */
static inline void
mac_to_crow(Macaulay* mac, const uint64_t idx, uint64_t* const mem,
            const uint64_t cr_slot_num, uint64_t* const drow,
            const uint64_t dr_slot_num) {
    memset(drow, 0x0, sizeof(uint64_t) * dr_slot_num);

    for(uint64_t i = 0; i < mac->eqs[idx].term_num; ++i) {
        drow_set_at(drow, mac->eqs[idx].row.s[i]);
    }

    uint64_t cursor = 0;
    for(uint64_t i = 0; i < dr_slot_num; ++i) {
        if(drow[i]) {
            *crow_idx(mem, cursor) = i;
            *crow_mono(mem, cursor) = drow[i];
            ++cursor;

            if(cursor >= cr_slot_num) {
                EXIT_WITH_MSG("[!] the ratio of present MQ terms is higher"
                              " than specified\n");
            }
        }
    }

    *crow_snum(mem) = cursor;

    mac->eqs[idx].term_num = MAC_CROW;
    mac->eqs[idx].row.c = mem;
}

/* function: mac_row_num
 * usage: given the number of variables and equations in the original MQ
 *      system, compute the number of rows in the Macaulay matrix.
 * arguments:
 *      1) eq_num: number of equations
 *      2) var_num: number of variables
 *      3) deg: degree of the Macaulay matrix
 * return: the number of rows
 */
inline uint64_t
mac_row_num(const uint64_t eq_num, const uint64_t var_num,
            const uint64_t deg) {
    assert(3 == deg || 4 == deg);
    uint64_t rol_num = (var_num + 1) * eq_num;

    if(4 == deg) {
        rol_num += binom(var_num, 2) * eq_num;
    }

    return rol_num;
}

/* function: mac_col_num
 * usage: given the number of variables and degree of the Macaulay matrix,
 *      compute the number of columns in the Macaulay matrix.
 * arguments:
 *      1) var_num: number of variables
 *      2) deg: degree of the Macaulay matrix
 * return: the number of columns
 */
inline uint64_t
mac_col_num(const uint64_t var_num, const uint64_t deg) {
    assert(3 == deg || 4 == deg);
    return sum_binom(var_num, deg);
}

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
mac_nl_num(const uint64_t var_num, const uint64_t deg, const uint64_t fix_var_num) {
    assert(3 == deg || 4 == deg);

    // quadratic terms
    uint64_t count = binom(var_num-fix_var_num, 2);

    // degree 3 monomials
    count += binom(var_num-fix_var_num, 3) +
        binom(var_num-fix_var_num, 2) * fix_var_num;

    if(4 == deg) {
        count += binom(var_num-fix_var_num, 4) +
            binom(var_num-fix_var_num, 3) * fix_var_num +
            binom(var_num-fix_var_num, 2) * binom(fix_var_num, 2);
    }

    return count;
}

/* A sub-routine of mac_init_perm, given square free equation, fill it into the
 * bottom of Macaulay matrix
 */
static inline void
cpy_eq_perm(MacRow* const __restrict__ eq, const bool* __restrict__ src,
            const uint64_t mq_term_num, const uint64_t mac_term_num,
            const uint64_t* const __restrict__ idx_map,
            uint64_t* const __restrict__ clz, bool* const __restrict__ flags) {

    memset(flags, 0x0, sizeof(bool) * mac_term_num);

    for(uint64_t i = 0; i < mq_term_num; ++i) {
        bool src_bit = src[mq_term_num-1-i];

        if(!src_bit) {
            continue;
        }
        
        uint64_t dest_idx = idx_map[mac_term_num-1-i];
        flags[dest_idx] = true;
    }

    for(uint64_t i = 0; i < mac_term_num; ++i) {
        if(flags[i]) {
            *clz = i;
            break;
        }
    }

    for(uint64_t i = *clz; i < mac_term_num; ++i) {
        if(!flags[i]) {
            continue;
        }
        
        eq->row.s[eq->term_num++] = i;
    }
}

/* function: cmp_indices
 * usage: compute the new indices of all monomial with degree <= 2 based on one
 *      multiplier
 * arguments:
 *      1) new_indices: container for the new indices
 *      2) var_num: number of variables in the system
 *      3) deg: degree of the Macaulay matrix
 *      4) var_idx: for multiplier x_{var_idx+1}
 * return: void
 */
static void
cmp_indices(uint64_t* new_indices, const uint64_t var_num, const uint64_t deg,
            const uint64_t var_idx) {
    
    uint64_t indices[4]; // at most 4 vars
    uint64_t src_idx = binom(var_num, 3);

    if(4 == deg) {
        src_idx += binom(var_num, 4);
    }

    // degree 2 monomial from src
    for(uint64_t i = 1; i < var_idx; ++i) {
        for(uint64_t j = 0; j < i; ++j) {
            indices[0] = j;
            indices[1] = i;
            indices[2] = var_idx;
            new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);
        }
    }

    for(uint64_t j = 0; j < var_idx; ++j) {
        new_indices[src_idx] = src_idx;
        ++src_idx;
    }

    for(uint64_t i = var_idx+1; i < var_num; ++i) {
        for(uint64_t j = 0; j < var_idx; ++j) {
            indices[0] = j;
            indices[1] = var_idx;
            indices[2] = i;
            new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);
        }

        new_indices[src_idx] = src_idx;
        ++src_idx;

        for(uint64_t j = var_idx+1; j < i; ++j) {
            indices[0] = var_idx;
            indices[1] = j;
            indices[2] = i;
            new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);
        }
    }

    // linear term from src
    for(uint64_t i = 0; i < var_idx; ++i) {
        indices[0] = i;
        indices[1] = var_idx;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 2);
    }

    new_indices[src_idx] = src_idx;
    ++src_idx;

    for(uint64_t i = var_idx+1; i < var_num; ++i) {
        indices[0] = var_idx;
        indices[1] = i;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 2);
    }

    // constant term from src
    indices[0] = var_idx;
    new_indices[src_idx++] = mono_idx(deg, var_num, indices, 1);
}


/* function: cmp_indices2
 * usage: compute the new indices of all monomial with degree <= 2 based on two
 *      multipliers
 * arguments:
 *      1) new_indices: container for the new indices
 *      2) var_num: number of variables in the system
 *      3) deg: degree of the Macaulay matrix
 *      4) min_var_idx: for multiplier x_{min_var_idx+1}, the smaller one
 *      5) max_var_idx: for multiplier x_{max_var_idx+1}, the larger one
 * return: void
 */
static void
cmp_indices2(uint64_t* new_indices, const uint64_t var_num, const uint64_t deg,
             const uint64_t min_var_idx, const uint64_t max_var_idx) {

    uint64_t src_idx = binom(var_num, 3);
    if(4 == deg) {
        src_idx += binom(var_num, 4);
    }

    uint64_t indices[4]; // at most degree 4
    
    // degree 2 monomials
    for(uint64_t i = 1; i < min_var_idx; ++i) {
        for(uint64_t j = 0; j < i; ++j) {
            indices[0] = j;
            indices[1] = i;
            indices[2] = min_var_idx;
            indices[3] = max_var_idx;
            new_indices[src_idx++] = mono_idx(deg, var_num, indices, 4);
        }
    }

    for(uint64_t j = 0; j < min_var_idx; ++j) {
        indices[0] = j;
        indices[1] = min_var_idx;
        indices[2] = max_var_idx;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);
    }

    for(uint64_t i = min_var_idx+1; i < max_var_idx; ++i) {
        for(uint64_t j = 0; j < min_var_idx; ++j) {
            indices[0] = j;
            indices[1] = min_var_idx;
            indices[2] = i;
            indices[3] = max_var_idx;
            new_indices[src_idx++] = mono_idx(deg, var_num, indices, 4);
        }

        indices[0] = min_var_idx;
        indices[1] = i;
        indices[2] = max_var_idx;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);

        for(uint64_t j = min_var_idx+1; j < i; ++j) {
            indices[0] = min_var_idx;
            indices[1] = j;
            indices[2] = i;
            indices[3] = max_var_idx;
            new_indices[src_idx++] = mono_idx(deg, var_num, indices, 4);
        }
    }

    // i = max_var_idx
    for(uint64_t j = 0; j < min_var_idx; ++j) {
        indices[0] = j;
        indices[1] = min_var_idx;
        indices[2] = max_var_idx;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);
    }

    // i = max_var_idx and j = min_var_idx
    new_indices[src_idx] = src_idx;
    ++src_idx;


    for(uint64_t j = min_var_idx+1; j < max_var_idx; ++j) {
        indices[0] = min_var_idx;
        indices[1] = j;
        indices[2] = max_var_idx;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);
    }

    for(uint64_t i = max_var_idx+1; i < var_num; ++i) {
        for(uint64_t j = 0; j < min_var_idx; ++j) {
            indices[0] = j;
            indices[1] = min_var_idx;
            indices[2] = max_var_idx;
            indices[3] = i;
            new_indices[src_idx++] = mono_idx(deg, var_num, indices, 4);
        }

        // j = min_var_idx
        indices[0] = min_var_idx;
        indices[1] = max_var_idx;
        indices[2] = i;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);

        for(uint64_t j = min_var_idx+1; j < max_var_idx; ++j) {
            indices[0] = min_var_idx;
            indices[1] = j;
            indices[2] = max_var_idx;
            indices[3] = i;
            new_indices[src_idx++] = mono_idx(deg, var_num, indices, 4);
        }

        // j = max_var_idx
        indices[0] = min_var_idx;
        indices[1] = max_var_idx;
        indices[2] = i;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);

        for(uint64_t j = max_var_idx+1; j < i; ++j) {
            indices[0] = min_var_idx;
            indices[1] = max_var_idx;
            indices[2] = j;
            indices[3] = i;
            new_indices[src_idx++] = mono_idx(deg, var_num, indices, 4);
        }
    }

    // linear terms
    for(uint64_t i = 0; i < min_var_idx; ++i) {
        indices[0] = i;
        indices[1] = min_var_idx;
        indices[2] = max_var_idx;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);
    }

    indices[0] = min_var_idx;
    indices[1] = max_var_idx;
    new_indices[src_idx++] = mono_idx(deg, var_num, indices, 2);

    for(uint64_t i = min_var_idx+1; i < max_var_idx; ++i) {
        indices[0] = min_var_idx;
        indices[1] = i;
        indices[2] = max_var_idx;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);
    }

    indices[0] = min_var_idx;
    indices[1] = max_var_idx;
    new_indices[src_idx++] = mono_idx(deg, var_num, indices, 2);

    for(uint64_t i = max_var_idx+1; i < var_num; ++i) {
        indices[0] = min_var_idx;
        indices[1] = max_var_idx;
        indices[2] = i;
        new_indices[src_idx++] = mono_idx(deg, var_num, indices, 3);
    }

    // constant
    indices[0] = min_var_idx;
    indices[1] = max_var_idx;
    new_indices[src_idx++] = mono_idx(deg, var_num, indices, 2);
}

/* A sub-routine for mac_init_perm */
static inline void
muleq_perm(MacRow* const __restrict__ dst, const bool* const __restrict__ src,
           const uint64_t mq_term_num, const uint64_t mac_term_num,
           const uint64_t* const __restrict__ indices,
           const uint64_t* const __restrict__ idx_map,
           uint64_t* const __restrict__ clz,
           bool* const __restrict__ flags) {

    memset(flags, 0x0, sizeof(bool) * mac_term_num);

    for(uint64_t i = 0; i < mq_term_num; ++i) {
        bool src_bit = src[mq_term_num-1-i];

        if(!src_bit) {
            continue;
        }
        
        uint64_t dest_idx = idx_map[indices[mac_term_num-1-i]];
        flags[dest_idx] = !flags[dest_idx];
    }

    for(uint64_t i = 0; i < mac_term_num; ++i) {
        if(flags[i]) {
            *clz = i;
            break;
        }
    }

    for(uint64_t i = *clz; i < mac_term_num; ++i) {
        if(!flags[i]) {
            continue;
        }
        
        dst->row.s[dst->term_num++] = i;
    }
}

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
              uint64_t* const __restrict__ clz) {
    assert(var_num == mac->var_num);
    assert(eq_num == mac->mq_eq_num);
    assert(3 == mac->deg || 4 == mac->deg);

    // init all entries to UINT64_MAX
    memset(clz, 0xFF, mac->eq_num * sizeof(uint64_t));
    bool* flags = SMALLOC(bool, mac->term_num);

    const uint64_t deg4_eq_offset = (3 == mac->deg) ? 0 : (binom(var_num, 2) * eq_num);

    // reduce the original system
    const uint64_t ori_eq_offset = deg4_eq_offset + var_num * eq_num;

    for(uint64_t eq_idx = 0; eq_idx < eq_num; ++eq_idx) {
        cpy_eq_perm(mac->eqs + ori_eq_offset + eq_idx, mqsys + term_num * eq_idx,
                    term_num, mac->term_num, idx_map, clz + ori_eq_offset + eq_idx,
                    flags);
    }

    uint64_t* indices = SMALLOC(uint64_t, mac->term_num);
    uint64_t dest_eq_idx = 0;

    if(4 == mac->deg) {
        // xixj * eq
        for(uint64_t var1_idx = 1; var1_idx < var_num; ++var1_idx) { // skip x1
            for(uint64_t var2_idx = 0; var2_idx < var1_idx; ++var2_idx) {

                cmp_indices2(indices, mac->var_num, mac->deg, var2_idx, var1_idx);

                // TODO: parallelize this?
                for(uint64_t eq_idx = 0; eq_idx < eq_num; ++eq_idx) {

                    muleq_perm(mac->eqs + dest_eq_idx, mqsys + term_num *  eq_idx,
                               term_num, mac->term_num, indices, idx_map,
                               clz + dest_eq_idx, flags);
                    ++dest_eq_idx;
                }
            }
        }
    }

    // xi * eq
    for(uint64_t var_idx = 0; var_idx < mac->var_num; ++var_idx) {
        cmp_indices(indices, mac->var_num, mac->deg, var_idx);

        // TODO: parallelize this?
        for(uint64_t eq_idx = 0; eq_idx < eq_num; ++eq_idx) {

            muleq_perm(mac->eqs + dest_eq_idx, mqsys + term_num * eq_idx,
                       term_num, mac->term_num, indices, idx_map,
                       clz + dest_eq_idx, flags);
            ++dest_eq_idx;
        }
    }
    
    free(indices);
    free(flags);
    assert(dest_eq_idx + eq_num == mac->eq_num);
}

/* function: mac_at
 * usage: return the entry given its position
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) row_idx: row index
 *      3) col_idx: column index
 * return: the bit at the given position
 */
inline bool
mac_at(Macaulay* mac, const uint64_t row_idx, const uint64_t col_idx) {

    if(MAC_DROW == mac->eqs[row_idx].term_num) {
        return drow_at(mac->eqs[row_idx].row.d, col_idx);
    }

    if(MAC_CROW == mac->eqs[row_idx].term_num) {
        uint64_t slot_idx = drow_slot_idx(col_idx);

        for(uint64_t i = 0; i < *crow_snum(mac->eqs[row_idx].row.c); ++i) {
            if( *crow_idx(mac->eqs[row_idx].row.c, i) == slot_idx ) {
                uint64_t mono = *crow_mono(mac->eqs[row_idx].row.c, i);
                return (mono >> drow_slot_offset(col_idx)) & 0x1UL;
            } else if( *crow_idx(mac->eqs[row_idx].row.c, i) > slot_idx) {
                return false;
            }
        }

        return false;
    }

    // TODO: optimize this with bsearch
    for(uint64_t i = 0; i < mac->eqs[row_idx].term_num; ++i) {
        if(col_idx == mac->eqs[row_idx].row.s[i]) {
            return true;
        } else if(col_idx < mac->eqs[row_idx].row.s[i]) {
            break;
        }
    }

    return false;
}

/* function: mac_zrow
 * usage: check if a row is all zeros
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) row_idx: row index
 * return: true if the row is all zero, otherwise false
 */
inline bool
mac_zrow(Macaulay* mac, const uint64_t row_idx) {
    if(MAC_CROW == mac->eqs[row_idx].term_num) {
        return (0 == *crow_snum(mac->eqs[row_idx].row.c));
    }

    if(MAC_DROW != mac->eqs[row_idx].term_num) {
        return (0 == mac->eqs[row_idx].term_num);
    }

    for(uint64_t i = 0; i < drow_slots_num(mac->term_num); ++i) {
        if( mac->eqs[row_idx].row.d[i] ) {
            return false;
        }
    }

    return true;
}

/* function: mac_row_clz
 * usage: count the number of leading zero monomials in a row
 * arguments:
 *      1) mac: pointer to structure Macaulay
 *      2) row_idx: row index
 * return: the number of leading zero monomials
 */
inline uint64_t
mac_row_clz(Macaulay* mac, const uint64_t row_idx) {
    if(MAC_CROW == mac->eqs[row_idx].term_num) {
        uint64_t lz = *crow_idx(mac->eqs[row_idx].row.c, 0) * 64;
        return __builtin_clzll(*crow_mono(mac->eqs[row_idx].row.c, 0)) + lz;
    }

    if(MAC_DROW != mac->eqs[row_idx].term_num) {
        return mac->eqs[row_idx].row.s[0];
    }

    uint64_t i;
    for(i = 0; i < drow_slots_num(mac->term_num); ++i) {
        if(mac->eqs[row_idx].row.d[i]) {
            break;
        }
    }

    return i * 64 + __builtin_clzll(mac->eqs[row_idx].row.d[i]);
}

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
mac_nl_indices(Macaulay* mac, uint64_t* indices, const uint64_t kvar_num) {
    uint64_t cursor = 0;
    uint64_t col_idx = 0;
    uint64_t bound = 0;

    if(4 == mac->deg) {  // degree 4 monomials
        for(uint64_t c = 0; c < binom(kvar_num, 4); ++c) {
            indices[cursor++] = col_idx++;
        }

        for(uint64_t i = kvar_num; i < mac->var_num; ++i) {
            for(uint64_t c = 0; c < binom(kvar_num, 3); ++c) {
                indices[cursor++] = col_idx++;
            }

            for(uint64_t j = kvar_num; j < i; ++j) {

                bound = binom(kvar_num, 2);
                for(uint64_t c = 0; c < bound; ++c) {
                    indices[cursor++] = col_idx++;
                }
             
                col_idx += binom(j, 2) - bound;
            }
        }
    }

    // degree 3 monomials
    for(uint64_t c = 0; c < binom(kvar_num, 3); ++c) {
        indices[cursor++] = col_idx++;
    }

    for(uint64_t i = kvar_num; i < mac->var_num; ++i) {
        bound = binom(kvar_num, 2);
        for(uint64_t c = 0; c < bound; ++c) {
            indices[cursor++] = col_idx++;
        }

        col_idx += binom(i, 2) - bound;
    }

    // degree 2 monomials
    bound = binom(kvar_num, 2);
    for(uint64_t c = 0; c < bound; ++c) {
        indices[cursor++] = col_idx++;
    }
}

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
                uint64_t* const __restrict__ l_indices, const uint64_t kvar_num) {
    uint64_t nl_cursor = 0;
    uint64_t l_cursor = 0;
    uint64_t col_idx = 0;

    if(4 == deg) {
        // degree 4 monomials
        for(uint64_t i = 3; i < kvar_num; ++i) {
            for(uint64_t j = 2; j < i; ++j) {
                for(uint64_t k = 1; k < j; ++k) {
                    for(uint64_t l = 0; l < k; ++l) {
                        nl_indices[nl_cursor++] = col_idx++;
                    }
                }
            }
        }

        for(uint64_t i = kvar_num; i < var_num; ++i) {
            for(uint64_t j = 2; j < kvar_num; ++j) {
                for(uint64_t k = 1; k < j; ++k) {
                    for(uint64_t l = 0; l < k; ++l) {
                        nl_indices[nl_cursor++] = col_idx++;
                    }
                }
            }

            for(uint64_t j = kvar_num; j < i; ++j) {
                for(uint64_t k = 1; k < kvar_num; ++k) {
                    for(uint64_t l = 0; l < k; ++l) {
                        nl_indices[nl_cursor++] = col_idx++;
                    }
                }

                for(uint64_t k = kvar_num; k < j; ++k) {
                    for(uint64_t l = 0; l < k; ++l) {
                        l_indices[l_cursor++] = col_idx++;
                    }
                }
            }
        }
    }
    
    // degree 3 monomials
    for(uint64_t i = 2; i < kvar_num; ++i) {
        for(uint64_t j = 1; j < i; ++j) {
            for(uint64_t k = 0; k < j; ++k) {
                nl_indices[nl_cursor++] = col_idx++;
            }
        }
    }

    for(uint64_t i = kvar_num; i < var_num; ++i) {
        for(uint64_t j = 1; j < kvar_num; ++j) {
            for(uint64_t k = 0; k < j; ++k) {
                nl_indices[nl_cursor++] = col_idx++;
            }
        }

        for(uint64_t j = kvar_num; j < i; ++j) {
            for(uint64_t k = 0; k < j; ++k) {
                l_indices[l_cursor++] = col_idx++;
            }
        }
    }

    // degree 2 monomials
    for(uint64_t i = 1; i < kvar_num; ++i) {
        for(uint64_t j = 0; j < i; ++j) {
            nl_indices[nl_cursor++] = col_idx++;
        }
    }

    for(uint64_t i = kvar_num; i < var_num; ++i) {
        for(uint64_t j = 0; j < i; ++j) {
            l_indices[l_cursor++] = col_idx++;
        }
    }

    // linear and constant terms
    for(uint64_t i = 0; i < var_num + 1; ++i) {
        l_indices[l_cursor++] = col_idx++;
    }
}

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
            const uint64_t nlt_num, const uint64_t lt_num) {

    for(uint64_t i = 0; i < nlt_num; ++i) {
        idx_map[nl_indices[i]] = i;
    }

    for(uint64_t i = 0; i < lt_num; ++i) {
        idx_map[l_indices[i]] = nlt_num + i;
    }
}

static inline void
swap_rows(Macaulay* mac, const uint64_t row1, const uint64_t row2) {
    MacRow tmp = mac->eqs[row1];
    mac->eqs[row1] = mac->eqs[row2];
    mac->eqs[row2] = tmp;
}

/* perform row reduction on a dense row with a sparse row */
static inline void
reduc_with_srow(Macaulay* mac, const uint64_t didx, const uint64_t sidx) {

    for(uint64_t i = 0; i < mac->eqs[sidx].term_num; ++i) {
        drow_toggle_at(mac->eqs[didx].row.d, mac->eqs[sidx].row.s[i]);
    }
}

/* perform row reduction on a dense row with a compact row */
static inline void
reduc_with_crow(Macaulay* mac, const uint64_t didx, const uint64_t cidx) {

    for(uint64_t i = 0; i < *crow_snum(mac->eqs[cidx].row.c); ++i) {
        uint64_t idx = *crow_idx(mac->eqs[cidx].row.c, i);
        uint64_t mono = *crow_mono(mac->eqs[cidx].row.c, i);
        mac->eqs[didx].row.d[idx] ^= mono;
    }
}

/* A sub-routine of mac_calc_rmac, the container(dst) must be all zero */
static inline void
pack_drow(uint64_t* const __restrict__ dst,
          const uint64_t* const __restrict__ drow,
          const uint32_t* const __restrict__ pidx,
          const uint64_t mpvt_num, const uint64_t nltnum, const uint64_t ltnum) {

    // missing pivots
    for(uint64_t i = 0; i < mpvt_num; ++i) {
        if(drow_at(drow, pidx[i])) {
            drow_set_at(dst, i);
        }
    }

    // linear terms
    for(uint64_t i = 0; i < ltnum; ++i) {
        if(drow_at(drow, nltnum + i)) {
            drow_set_at(dst, mpvt_num + i);
        }
    }
}

/* arguments for worker thread of mac_calc_rmac */
typedef  struct {
    Macaulay* mac;
    const uint32_t* drows;
    uint64_t start;
    uint64_t piv;
    uint64_t begin;
    uint64_t end;
} crmac_arg;

/* worker thread of mac_calc_rmac */
static void
reduc_drow(void* dummy) {
    crmac_arg* arg = (crmac_arg*) dummy;

    for(uint64_t i = arg->begin; i < arg->end; ++i) {
        uint64_t eq = arg->drows[i];

        for(uint64_t j = arg->start; j < arg->piv; ++j) {
            if(drow_at(arg->mac->eqs[eq].row.d, j)) {
                reduc_with_crow(arg->mac, eq, j);
            }
        }
    }
}

/* arguments for worker thread of mac_calc_rmac */
typedef struct {
    Macaulay* mac;
    const uint32_t* drows;
    uint64_t drow_num;
    RMac* rmac;
    uint64_t ltnum;
    uint64_t begin;
    uint64_t end;
} cdr_arg;

/* worker thread of mac_calc_rmac */
static void
compact_drow(void* dummy) {
    cdr_arg* arg = (cdr_arg*) dummy;

    for(uint64_t i = arg->begin; i < arg->end; ++i) {
        pack_drow(rmac_row(arg->rmac, i), arg->mac->eqs[arg->drows[i]].row.d,
                  arg->drows, arg->drow_num, arg->mac->term_num - arg->ltnum,
                  arg->ltnum);
    }
}

/* arguments for worker thread of mac_calc_rmac */
typedef struct {
    Macaulay* mac;
    const uint32_t* drows;
    uint64_t start;
    uint64_t end;
    uint64_t nltnum;
} alf_arg;

/* worker thread of mac_calc_rmac */
static void
reduc_drow_all_found(void* dummy) {
    alf_arg* arg = (alf_arg*) dummy;

    for(uint64_t d = arg->start; d < arg->end; ++d) {
        uint64_t eq = arg->drows[d];

        // check each non-linear term in the eq
        for(uint64_t i = 0; i < arg->nltnum; ++i) {
            if(drow_at(arg->mac->eqs[eq].row.d, i)) {
                reduc_with_crow(arg->mac, eq, i);
            }
        }
    }
}

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
              threadpool_t* const tpool) {
    
    crmac_arg* const argpool = SMALLOC(crmac_arg, tpool->thread_count);
    cdr_arg* const cargpool = SMALLOC(cdr_arg, tpool->thread_count);
    alf_arg* afargpool = SMALLOC(alf_arg, tpool->thread_count);
    uint64_t chunk_size = (drow_num + keq_num) /  tpool->thread_count;

    if(drow_num) {
        // eliminate pivot terms that were found during pivoting
        uint64_t start = 0;
        for(uint64_t d = 0; d < drow_num; ++d) {
            uint64_t piv = drows[d];
            
            for(int i = 0; i < tpool->thread_count-1; ++i) {
                argpool[i].mac = mac;
                argpool[i].drows = drows;
                argpool[i].start = start;
                argpool[i].piv = piv;
                argpool[i].begin = i * chunk_size;
                argpool[i].end = (i+1) * chunk_size;

                if(threadpool_add(tpool, reduc_drow, argpool + i, 0)) {
                    EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
                }
            }

            argpool[tpool->thread_count-1].mac = mac;
            argpool[tpool->thread_count-1].drows = drows;
            argpool[tpool->thread_count-1].start = start;
            argpool[tpool->thread_count-1].piv = piv;
            argpool[tpool->thread_count-1].begin = (tpool->thread_count-1) * chunk_size;
            argpool[tpool->thread_count-1].end = drow_num + keq_num;

            if(threadpool_add(tpool, reduc_drow, argpool + (tpool->thread_count-1), 0)) {
                EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
            }

            if(threadpool_join(tpool, 0)) {
                EXIT_WITH_MSG("[!] failed to join thread pool\n");
            }

            start = piv + 1;
        }
    } else { // all pivot rows were found
        chunk_size = s->eq_num / tpool->thread_count;

        for(int i = 0; i < tpool->thread_count-1; ++i) {
            afargpool[i].mac = mac;
            afargpool[i].drows = drows;
            afargpool[i].start = i * chunk_size;
            afargpool[i].end = (i+1) * chunk_size;
            afargpool[i].nltnum = mac->term_num - ltnum;

            if(threadpool_add(tpool, reduc_drow_all_found, afargpool + i, 0)) {
                EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
            }
        }

        afargpool[tpool->thread_count-1].mac = mac;
        afargpool[tpool->thread_count-1].drows = drows;
        afargpool[tpool->thread_count-1].start = (tpool->thread_count-1) * chunk_size;
        afargpool[tpool->thread_count-1].end = s->eq_num;
        afargpool[tpool->thread_count-1].nltnum = mac->term_num - ltnum;

        if(threadpool_add(tpool, reduc_drow_all_found,
                          afargpool + tpool->thread_count-1, 0)) {
            EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
        }

        if(threadpool_join(tpool, 0)) {
            EXIT_WITH_MSG("[!] failed to join thread pool\n");
        }
    }

    chunk_size = (drow_num + keq_num) /  tpool->thread_count;
    // pact dense rows into a matrix of smaller dimensions
    for(int i = 0; i < tpool->thread_count-1; ++i) {
        cargpool[i].mac = mac;
        cargpool[i].drows = drows;
        cargpool[i].drow_num = drow_num;
        cargpool[i].rmac = s;
        cargpool[i].ltnum = ltnum;
        cargpool[i].begin = i * chunk_size;
        cargpool[i].end = (i+1) * chunk_size;

        if(threadpool_add(tpool, compact_drow, cargpool + i, 0)) {
            EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
        }
    }

    cargpool[tpool->thread_count-1].mac = mac;
    cargpool[tpool->thread_count-1].drows = drows;
    cargpool[tpool->thread_count-1].drow_num = drow_num;
    cargpool[tpool->thread_count-1].rmac = s;
    cargpool[tpool->thread_count-1].ltnum = ltnum;
    cargpool[tpool->thread_count-1].begin = (tpool->thread_count-1) * chunk_size;
    cargpool[tpool->thread_count-1].end = drow_num + keq_num;

    if(threadpool_add(tpool, compact_drow, cargpool + (tpool->thread_count-1), 0)) {
        EXIT_WITH_MSG("[!] failed to add worker to thread pool\n");
    }

    if(threadpool_join(tpool, 0)) {
        EXIT_WITH_MSG("[!] failed to join thread pool\n");
    }

    free(argpool);
    free(cargpool);
    free(afargpool);
}

/* function: mac_lzmap
 * usage: given the number of leading zeros of each row, compute the map from
 *      the number of leading zeros to the row index with that property. If no
 *      such row is found, UINT64_MAX will be filled instead. The domain is 0 ~ (n-1)
 * arguments:
 *      1) lzmap: container for the indices
 *      2) n: the domain of the map to compute, also the size of lzmap
 *      3) clz: the indices of the first non-zero entry of each row
 *      4) eq_num: the number of eqs, also the size of clz
 * return: the number of values in the domain without a suitable row
 */
static inline void
mac_lzmap(uint64_t* const __restrict__ lzmap, const uint64_t n,
          const uint64_t* const __restrict__ clz, const uint64_t eq_num) {
    memset(lzmap, 0xFF, n * sizeof(uint64_t));

    for(uint64_t i = 0; i < eq_num; ++i) {
        uint64_t lznum = clz[i];
        if(lznum < n) {
            lzmap[lznum] = i; // overwrite is fine
        }
    }
}

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
        uint32_t* const __restrict__ drows, const double mq_term_ratio) {

    uint64_t* lzmap = SMALLOC(uint64_t, col_num);
    uint64_t* crows = SMALLOC(uint64_t, col_num);
    mac_lzmap(lzmap, col_num, clz, mac->eq_num);

    uint64_t not_found = 0;
    uint64_t found = 0;
    for(uint64_t i = 0; i < col_num; ++i) {
        
        // find the perfect row
        uint64_t src_ridx = lzmap[i];
        if(UINT64_MAX == src_ridx) {
            drows[not_found++] = i;
            continue;
        }

        swap_rows(mac, i, src_ridx);
        
        // update lzmap and clz
        lzmap[i] = i;

        if(clz[i] < col_num) {
            lzmap[clz[i]] = src_ridx;
        }

        clz[src_ridx] = clz[i];
        clz[i] = i;

        crows[found++] = i;
    }

    const uint64_t dr_slot_num = drow_slots_num(mac->term_num);
    mac->dmem_blk = (uint64_t*) calloc((not_found + keq_num) * dr_slot_num,
                                       sizeof(uint64_t));
    if(NULL == mac->dmem_blk) {
        EXIT_WITH_MSG("[!] insufficient memory\n");
    }
    
    // transform rows that are not in their final position into dense format
    for(uint64_t i = 0; i < not_found; ++i) {
        mac_to_drow(mac, drows[i], mac->dmem_blk + dr_slot_num * i);
    }

    // randomly choose the equations to keep
    const uint64_t range = mac->eq_num - col_num;
    for(uint64_t i = 0; i < keq_num; ++i) {
        uint64_t idx = col_num + (rand() % range);
        swap_rows(mac, col_num + i, idx);
    }

    // transform the sub-system candidates into dense format
    for(uint64_t i = 0; i < keq_num; ++i) {
        mac_to_drow(mac, col_num + i,
                    mac->dmem_blk + dr_slot_num * (not_found + i));
        drows[not_found + i] = col_num + i;
    }

    // transform rows that are in their final position into compact format
    const uint64_t mq_tnum = binom(mac->var_num, 2) + mac->var_num + 1;
    const uint64_t cr_slot_num = crow_slot_num(mq_tnum, mq_term_ratio);
    mac->cmem_blk = (uint64_t*) malloc(mac_crow_memsize(mq_tnum, found, mq_term_ratio));
    if(NULL == mac->cmem_blk) {
        EXIT_WITH_MSG("[!] insufficient memory\n");
    }

    uint64_t* tmp = SMALLOC(uint64_t, dr_slot_num);
    for(uint64_t i = 0; i < found; ++i) {
        mac_to_crow(mac, crows[i], mac->cmem_blk + cr_slot_num * i,
                    cr_slot_num, tmp, dr_slot_num);
    }

    free(mac->mem_blk);
    mac->mem_blk = NULL;

    free(tmp);
    free(lzmap);
    free(crows);
    return not_found;
}

/* function: mono_idx
 * usage: given the degree and the variables in a monomial, computes its index
 *      according to grlex order.
 * arguments:
 *      1) deg: degree of the Macaulay matrix
 *      2) var_num: number of variables in the system
 *      3) vars: a sorted array of integers representing the monomial,
 *              e.g. for x1x3x4, one must pass [0, 2, 3]
 *      4) mvar_num: number of variables in the monomial
 * return: the index
 */
inline uint64_t
mono_idx(const uint64_t deg, const uint64_t var_num,
         const uint64_t* vars, const uint64_t mvar_num) {
    assert(mvar_num <= var_num);

    uint64_t idx = -1; // error, become extremely large uint64_t
    switch(mvar_num) {
        case 0: // constant term
            idx = binom(var_num, 4) + binom(var_num, 3) + binom(var_num, 2) + var_num;
            break;
        case 1: // linear term
            idx = binom(var_num, 4) + binom(var_num, 3) + binom(var_num, 2) + vars[0];
            break;
        case 2: // quadratic term
            idx = binom(var_num, 4) + binom(var_num, 3) + binom(vars[1], 2) + vars[0];
            break;
        case 3: // degree 3 monomial
            idx = binom(var_num, 4) + binom(vars[2], 3) + binom(vars[1], 2) + vars[0];
            break;
        case 4: // degree 4 monomial
            idx = binom(vars[3], 4) + binom(vars[2], 3) + binom(vars[1], 2) + vars[0];
            break;
        default:
            // do nothing
            break;
    }

    if(3 == deg) {
        idx -= binom(var_num, 4);
    }

    return idx;
}

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
inline uint64_t
fmono_idx(const uint64_t deg, const uint64_t var_num,
          const uint64_t* const vars, const uint64_t mvar_num) {
    assert(1 == mvar_num || 2 == mvar_num);

    uint64_t indices[4]; // at most degree 4
    indices[0] = 0; // x1
    indices[1] = 1; // x2
    indices[2] = 2; // x3
    indices[3] = 3; // x4
    uint64_t mdeg = 3;

    // the smaller multiplier
    if(2 < vars[0]) { // not x1 ~ x3
        indices[2] = vars[0];
    }

    // the larger multiplier
    if(2 == mvar_num) {
        if(3 < vars[1]) { // not x1 ~ x4
            indices[3] = vars[1];
        }
        ++mdeg;
    }
    
    return mono_idx(deg, var_num, indices, mdeg);
}

/* subroutine of mac_vbound */
// TODO: optimize this
static inline void
idx2vars(uint64_t* const vars, const uint64_t idx) {
    uint64_t i = 1;
    uint64_t offset = 0;

    while( (offset + i) < idx) {
        offset += i;
        ++i;
    }

    vars[1] = i;
    vars[0] = idx - offset;
}

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
inline uint64_t
mac_vbound(const uint64_t eq_num, const uint64_t var_num, const uint64_t deg,
           const uint64_t idx) {
    uint64_t vars[4];
    uint64_t eqgrp_idx = idx / eq_num;

    if(3 == deg) {
        if(eqgrp_idx < var_num) {
            vars[0] = eqgrp_idx;
            return fmono_idx(3, var_num, vars, 1);
        }
        
        // original MQ system
        vars[0] = 0;
        vars[1] = 1;
        return mono_idx(3, var_num, vars, 2); // x1x2
    }
    
    const uint64_t d2_sep = binom(var_num, 2);
    const uint64_t d1_sep = d2_sep + var_num;

    if(eqgrp_idx < d2_sep) {

        // xixj * MQ system
        idx2vars(vars, eqgrp_idx);
        return fmono_idx(4, var_num, vars, 2);

    } else if(eqgrp_idx < d1_sep) {

        // xi * MQ system
        vars[0] = eqgrp_idx - d2_sep;
        return fmono_idx(4, var_num, vars, 1);

    } else if(eqgrp_idx == d1_sep) {

        // original MQ system
        vars[0] = 0;
        vars[1] = 1;
        return mono_idx(4, var_num, vars, 2); // x1x2

    } else {
        return -1; // sth wrong, let it explode
    }
}

/* function: mac_hsnum
 * usage: compute the number of column where the horizontal boundary changes
 * arguments:
 *      1) var_num: number of variables
 *      2) deg: degree of the Macaulay matrix
 * return: the number
 */
inline uint64_t
mac_hsnum(const uint64_t var_num, const uint64_t deg) {
    uint64_t num = var_num - 1;

    if(4 == deg) {
        num += binom(var_num-2, 2);
    }

    return num;
}

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
          const uint64_t eq_num, const uint64_t var_num, const uint64_t deg) {
    uint64_t pass = 0;
    uint64_t vars[4];

    if(4 == deg) {
        ses[pass] = 0; // x1x2x3x4
        hbs[pass++] = 6 * eq_num; // choose(4, 2)

        for(uint64_t i = 4; i < var_num; ++i) {
            vars[0] = 0;
            vars[1] = 1;
            vars[2] = 2;
            vars[3] = i;
            ses[pass] = mono_idx(4, var_num, vars, 4); // x1x2x3x?
            hbs[pass] = hbs[pass-1] + 3 * eq_num;
            ++pass;

            for(uint64_t j = 3; j < i; ++j) {
                vars[2] = j;
                ses[pass] = mono_idx(4, var_num, vars, 4); // x1x2x?x?
                hbs[pass] = hbs[pass-1] + eq_num;
                ++pass;
            }
        }
    }

    vars[0] = 0;
    vars[1] = 1;
    vars[2] = 2;
    ses[pass] = mono_idx(deg, var_num, vars, 3); // x1x2x3
    hbs[pass] = 3 * eq_num + ( (3 == deg) ? 0 : hbs[pass-1] );
    ++pass;

    for(uint64_t i = 3; i < var_num; ++i) {
        vars[2] = i;
        ses[pass] = mono_idx(deg, var_num, vars, 3); // x1x2x?
        hbs[pass] = hbs[pass-1] + eq_num;
        ++pass;
    }

    ses[pass] = mono_idx(deg, var_num, vars, 2); // x1x2
    hbs[pass] = hbs[pass-1] + eq_num;
}


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
           const uint64_t* const __restrict__ hbs) {
    // TODO: use binary search
    uint64_t i = 0;
    while(i < senum && ses[i] < idx) {
        ++i;
    }

    return hbs[i];
}
