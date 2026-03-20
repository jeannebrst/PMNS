#ifndef MULT_15_H
#define MULT_15_H

#include <stdint.h>
#include <inttypes.h>
#include <arm_neon.h>
#include "mphell/mphell.h"

static amns_block gen_Mti15[29];
static amns_block gen_Mtiprime15[29];
static amns_block gen_B[29];

void get_gen_from_vec(amns_block gen[29], amns_block B[15], amns_block lambda);
void get_gen_15(amns_block gen_Mti[29], amns_elt * Mti);

void mult_jeanne(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[29], amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]);
void mult_jeanne_neon(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[29], amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]);

#endif