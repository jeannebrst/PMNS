#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <arm_neon.h>

#include "mult_16.h"
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

    // char *amns_param = get_amns_param(char_p, "16");
    char *amns_param = "[15, 16, [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 18, 1732996756754789471339355312096369946981361386195029723648, [-2239, -494, -130, -54, -772, 615, -474, -1531, 1849, -895, -1125, -289, 438, -1978, 590, 1077], [4092435061, 3150241284, 3437500926, 2483521910, 1913883334, 2800600421, 2868738738, 3215995195, 1101386429, 229037773, 2613946518, 4042815467, 3393698058, 1925116940, 1693409253, 1245900174], [-5701, -3453, -5708, -4336, -4736, -3772, -3384, -5334, -4480, -2794, -2698, -4467, -761, -3023, -3378, -1571], [-1812, -10495, 232, -4276, -8242, -1249, -1563, -6222, -3852, 97, -6506, -2222, 120, -4246, -3666, -1254]]";

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

    get_gen_16(gen_Mti16, ((fp_param)((field_srcptr)f)->param)->AMNS->Mti);
    get_gen_16(gen_Mtiprime16, ((fp_param)((field_srcptr)f)->param)->AMNS->Mtiprime);
    get_gen_from_vec(gen_B, *(fp_elt_srcptr)(*b), ((fp_param)((field_srcptr)f)->param)->AMNS->lambda);

    int rec = 2; // paramètre pour niveau de récursivité de certaines opérations, peut être 0, 1 ou 2 (si > 2, sera interprété comme 2)

    mult_jeanne(c, a, b, f, STACK_1, gen_B, gen_Mtiprime16, gen_Mti16, rec);
    printf("\nNormal : a*b = "); field_elt_print(c, 10, true, f, STACK_1); printf("\n\n");

    mult_jeanne_neon(c, a, b, f, STACK_1, gen_B, gen_Mtiprime16, gen_Mti16, rec);
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