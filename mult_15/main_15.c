#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <arm_neon.h>

#include "mult_15.h"
#include "../util.h"
#include "mphell/mphell.h"

int main(int argc, char **argv) {

    char char_a[1024], char_b[1024], char_p[1024];

    if (argc == 4) {
        snprintf(char_a, sizeof(char_a), "%s", argv[1]);
        snprintf(char_b, sizeof(char_b), "%s", argv[2]);
        snprintf(char_p, sizeof(char_p), "%s", argv[3]);
    } else {
        char input[1024];

        printf("a b p ? (or 'd' for default values)\n");
        scanf("%1023s", input);

        if (strcmp(input, "d") == 0) {
            snprintf(char_a, sizeof(char_a), "10");
            snprintf(char_b, sizeof(char_b), "10");
            snprintf(char_p, sizeof(char_p), "6277101735386680763835789423207666416083908700390324961279");
        } else {
            snprintf(char_a, sizeof(char_a), "%s", input);
            scanf("%1023s %1023s", char_b, char_p);
        }
    }

    init_mphell(256, RANDOM_AES256, DEVURANDOM);

    field f;
    field_alloc(f, FP, bits_to_nblock(256), NULL);

    number p;
    number_init(&p, bits_to_nblock(256));

    number_set_str(p, char_p, 10);

    char *amns_param = get_amns_param(char_p, "15");
    // char *amns_param = "[15, 15, [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 18, 4647938237712258623232008979259790892265608708830902830182, [-1715, 2144, 90, -1538, 2908, -3091, -277, 701, -150, -1017, 2536, 2996, -2320, 191, -3131], [1419663213, 1657114742, 3852025598, 3602762570, 3511206674, 3088878177, 161064579, 3194799017, 1912074464, 2361726565, 1147698033, 1518907444, 166455973, 301739933, 1505293996], [-5334, 3887, -7054, -6069, 2955, -1940, -4097, 8067, -2335, -5784, 3054, -6574, -636, 1487, -2157], [-4194, -4814, -2611, -2853, -3020, 405, 1338, 1791, -2395, -4058, -1256, -3941, -1570, -2012, 2154]]";
    
    amns AMNS;

    amns_alloc_init_str(&AMNS, amns_param, p);

    field_set_amns(f, AMNS);
    field_create(f, "", STACK_1, 1, p);

    field_elt a, b, c;
    field_elt_alloc(&a, f);
    field_elt_init(a, f);
    field_elt_alloc(&b, f);
    field_elt_init(b, f);
    field_elt_alloc(&c, f);
    field_elt_init(c, f);

    field_elt_set_str(a, char_a, 10, false, f, STACK_1);
    field_elt_set_str(b, char_b, 10, false, f, STACK_1);

    printf("a = "); field_elt_print(a, 10, true, f, STACK_1); printf("\n");
    printf("b = "); field_elt_print(b, 10, true, f, STACK_1); printf("\n");

    get_gen_15(gen_Mti15, ((fp_param)((field_srcptr)f)->param)->AMNS->Mti);
    get_gen_15(gen_Mtiprime15, ((fp_param)((field_srcptr)f)->param)->AMNS->Mtiprime);
    get_gen_from_vec(gen_B, *(fp_elt_srcptr)(*b), ((fp_param)((field_srcptr)f)->param)->AMNS->lambda);

    mult_jeanne(c, a, b, f, STACK_1, gen_B, gen_Mtiprime15, gen_Mti15);
    printf("\nNormal : a*b = "); field_elt_print(c, 10, true, f, STACK_1); printf("\n\n");

    mult_jeanne_neon(c, a, b, f, STACK_1, gen_B, gen_Mtiprime15, gen_Mti15);
    printf("\nNEON : a*b = "); field_elt_print(c, 10, true, f, STACK_1); printf("\n\n");

    if(amns_param) {
        printf("Sortie de sage :\n%s", amns_param);
        printf("\n");
        //free(amns_param);
    } else {
        printf("Erreur lors de l'exécution sage\n");
    }

    return 0;
}