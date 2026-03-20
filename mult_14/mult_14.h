#ifndef MULT_14_H
#define MULT_14_H

#include <stdint.h>
#include <inttypes.h>
#include <arm_neon.h>
#include <mphell/mphell.h>
#include <mphell/mphell-amns.h>

static amns_block gen_Mti14[27];
static amns_block gen_Mtiprime14[27];
static amns_block gen_B[27];

void get_gen_from_vec(amns_block gen[27], amns_block B[14], amns_block lambda);
void get_gen_14(amns_block gen_Mti[27], amns_elt * Mti);

void mult_jeanne(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[27], amns_block gen_Mtiprime14[27], amns_block gen_Mti14[27]);
void mult_jeanne_neon(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[27], amns_block gen_Mtiprime14[27], amns_block gen_Mti14[27]);

#endif