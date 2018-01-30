/* drow.h: header file of functions for manipulating drows.
 */

#ifndef __DROW_H__
#define __DROW_H__

#include <stdint.h>
#include <stdbool.h>


/* function: drow_slots_num
 * usage: compute how many slots are needed
 * argument: mnum, number of monomials
 * return: the number of slots
 */
uint64_t
drow_slots_num(const uint64_t mnum);

/* function: drow_slot_idx
 * usage: given the monomial index, compute the index of the slot where it locates
 * arguments:
 *      1) idx: index of the monomial
 * return: the slot index
 */
uint64_t 
drow_slot_idx(const uint64_t idx);

/* function: drow_slot_offset
 * usage: given the monomial index, compute the offset from RHS into slot where
 *      it belongs
 * arguments:
 *      1) i: index of the monomial
 * return: the slot index
 */
uint64_t
drow_slot_offset(const uint64_t i);

/* function: drow_set_at
 * usage: set the monomial at the given index to 1
 * arguments:
 *      1) drow: a pointer to the dense row
 *      2) i: index of the monomial, 0 is first monomial is grlex
 * return: void
 */
void
drow_set_at(uint64_t* const drow, uint64_t i);

/* function: drow_at
 * usage: return the monomial at the given index
 * arguments:
 *      1) drow: a pointer to the dense row
 *      2) i: index of the monomial, 0 is first monomial is grlex
 * return: void
 */
bool
drow_at(const uint64_t* const drow, uint64_t i);

/* function: drow_toggle_at
 * usage: toggle the monomial at the given index
 * arguments:
 *      1) drow: a pointer to the dense row
 *      2) i: index of the monomial, 0 is first monomial is grlex
 * return: void
 */
void
drow_toggle_at(uint64_t* const drow, uint64_t i);

#endif
