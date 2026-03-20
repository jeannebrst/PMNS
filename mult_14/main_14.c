#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <arm_neon.h>

#include "mult_14.h"
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
            snprintf(char_p, sizeof(char_p), "115792089210356248762697446949407573530086143415290314195533631308867097853951");
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

    char *amns_param = get_amns_param(char_p, "14");
    // char *amns_param = "[1, 14, [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 24, 114138283086017516355915752550460811928665718035587855046549153961407115672824, [4683, 48932, -188906, -57287, -6438, 137188, -128445, -71221, 14846, 154270, -153728, -45026, 105602, -68045], [4193936703, 2542287792, 4133871222, 1756401359, 1193602640, 1725648838, 2893131054, 3894537603, 3759292924, 2109451957, 2628139664, 2676699888, 1222600349, 1316339254], [-203232, 84149, -397332, -288348, -125646, -103816, -20848, -480462, -99047, -65942, -113492, -80851, -262643, -230489], [-140679, -148002, -156171, -482730, 120772, -44807, -511768, -256266, 170912, -323770, -145291, -130161, -67380, -155425]]";

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

    get_gen_14(gen_Mti14, ((fp_param)((field_srcptr)f)->param)->AMNS->Mti);
    get_gen_14(gen_Mtiprime14, ((fp_param)((field_srcptr)f)->param)->AMNS->Mtiprime);
    get_gen_from_vec(gen_B, *(fp_elt_srcptr)(*b), ((fp_param)((field_srcptr)f)->param)->AMNS->lambda);

    mult_jeanne(c, a, b, f, STACK_1, gen_B, gen_Mtiprime14, gen_Mti14);
    printf("\nNormal : a*b = "); field_elt_print(c, 10, true, f, STACK_1); printf("\n\n");

    mult_jeanne_neon(c, a, b, f, STACK_1, gen_B, gen_Mtiprime14, gen_Mti14);
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