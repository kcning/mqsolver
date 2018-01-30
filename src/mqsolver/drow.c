/* drow.c: implementation for header file drow.h */

#include "drow.h"

#include <stdint.h>
#include <stdbool.h>


/* function: drow_slots_num
 * usage: compute how many slots are needed
 * argument: mnum, number of monomials
 * return: the number of slots
 */
inline uint64_t
drow_slots_num(const uint64_t mnum) {
    return (1 + ((mnum - 1) / 64));
}

/* function: drow_slot_idx
 * usage: given the monomial index, compute the index of the slot where it locates
 * arguments:
 *      1) idx: index of the monomial
 * return: the slot index
 */
inline uint64_t 
drow_slot_idx(const uint64_t idx) {
    return idx >> 6;  /* divide by 64 */
}

/* function: drow_slot_offset
 * usage: given the monomial index, compute the offset from RHS into slot where
 *      it belongs
 * arguments:
 *      1) i: index of the monomial
 * return: the slot index
 */
inline uint64_t
drow_slot_offset(const uint64_t i) {
    return 63 - (i & 0x3FUL); // module 64, then reverse
}

/* function: drow_set_at
 * usage: set the monomial at the given index to 1
 * arguments:
 *      1) drow: a pointer to the dense row
 *      2) i: index of the monomial, 0 is first monomial is grlex
 * return: void
 */
inline void
drow_set_at(uint64_t* const drow, uint64_t i) {
    drow[drow_slot_idx(i)] |= 0x1UL << drow_slot_offset(i);
}

/* function: drow_at
 * usage: return the monomial at the given index
 * arguments:
 *      1) drow: a pointer to the dense row
 *      2) i: index of the monomial, 0 is first monomial is grlex
 * return: void
 */
inline bool
drow_at(const uint64_t* const drow, uint64_t i) {
    return (drow[drow_slot_idx(i)] >> drow_slot_offset(i)) & 0x1UL;
}

/* function: drow_toggle_at
 * usage: toggle the monomial at the given index
 * arguments:
 *      1) drow: a pointer to the dense row
 *      2) i: index of the monomial, 0 is first monomial is grlex
 * return: void
 */
inline void
drow_toggle_at(uint64_t* const drow, uint64_t i) {
    drow[drow_slot_idx(i)] ^= 0x1UL << drow_slot_offset(i);
}
