#ifndef MULT_16_H
#define MULT_16_H

#include <stdint.h>
#include <inttypes.h>
#include <arm_neon.h>
#include "mphell/mphell.h"

static amns_block gen_Mti16[31];
static amns_block gen_Mtiprime16[31];
static amns_block gen_B[31];

void get_gen_from_vec(amns_block gen[31], amns_block B[16], amns_block lambda);
void get_gen_16(amns_block gen_Mti[31], amns_elt * Mti);

void mult_jeanne(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[31], amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]);
void mult_jeanne_neon(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[31], amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]);

#endif