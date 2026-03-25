#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <arm_neon.h>

#include "mult_16.h"
#include "../util.h"
#include "mphell/mphell.h"

#define BILLION 1000000000L

int main(int argc, char **argv) {

    init_mphell(256, RANDOM_AES256, DEVURANDOM);
    int nb_iteration = 1<<20;
    printf("Number of iteration: %d\n", nb_iteration); 
    int i;
    bool test = true;
    struct timespec start1, end1;
    long unsigned int diff1, diff2, diff3;

    field f;
    field_alloc(f, FP, bits_to_nblock(256), NULL);

    number p;
    number_init(&p, bits_to_nblock(256));
    number_set_str(p, "6277101735386680763835789423207666416083908700390324961279", 10);

    // char *amns_param = get_amns_param("6277101735386680763835789423207666416083908700390324961279", "16");
    char *amns_param = "[15, 16, [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 18, 1732996756754789471339355312096369946981361386195029723648, [-2239, -494, -130, -54, -772, 615, -474, -1531, 1849, -895, -1125, -289, 438, -1978, 590, 1077], [4092435061, 3150241284, 3437500926, 2483521910, 1913883334, 2800600421, 2868738738, 3215995195, 1101386429, 229037773, 2613946518, 4042815467, 3393698058, 1925116940, 1693409253, 1245900174], [-5701, -3453, -5708, -4336, -4736, -3772, -3384, -5334, -4480, -2794, -2698, -4467, -761, -3023, -3378, -1571], [-1812, -10495, 232, -4276, -8242, -1249, -1563, -6222, -3852, 97, -6506, -2222, 120, -4246, -3666, -1254]]";
    
    amns AMNS;

    amns_alloc_init_str(&AMNS, amns_param, p);

    // printf("The AMNS system is:\n"); amns_print_AMNS(AMNS); printf("\n");

    field_set_amns(f, AMNS);
    field_create(f, "", STACK_1, 1, p);

    number a;
    number b;
    field_elt c;
    field_elt d;
    field_elt res;

    number_init(&a, bits_to_nblock(AMNS->p_size));
    number_init(&b, bits_to_nblock(AMNS->p_size));

    number_random1(a, p, STACK_1);
    number_random1(b, p, STACK_1);

    printf("a = "); number_print(a, 10); printf("\n");
    printf("b = "); number_print(b, 10); printf("\n");

    field_elt_alloc(&c, f);
    field_elt_init(c, f);
    field_elt_alloc(&d, f);
    field_elt_init(d, f);
    field_elt_alloc(&res, f);
    field_elt_init(res, f);

    field_elt_set_number(c, false, f, STACK_1, 1, a);
    field_elt_set_number(d, false, f, STACK_1, 1, b);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
    for (i=0; i<nb_iteration; i++)
    {
        field_elt_mul(res, c, d, f, STACK_1);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
    diff1 = BILLION * (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec);

    printf("mean time using field_elt_mul() = \t%lu nanoseconds\n", diff1 / nb_iteration);

    get_gen_16(gen_Mti16, ((fp_param)((field_srcptr)f)->param)->AMNS->Mti);
    get_gen_16(gen_Mtiprime16, ((fp_param)((field_srcptr)f)->param)->AMNS->Mtiprime);
    get_gen_from_vec(gen_B, *(fp_elt_srcptr)(*d), ((fp_param)((field_srcptr)f)->param)->AMNS->lambda);

    int rec = 2; // paramètre pour niveau de récursivité de certaines opérations, peut être 0, 1 ou 2 (si > 2, sera interprété comme 2)

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
    for (i=0; i<nb_iteration; i++)
    {
        mult_jeanne(res, c, d, f, STACK_1, gen_B, gen_Mtiprime16, gen_Mti16, rec);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
    diff2 = BILLION * (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec);

    printf("mean time using mult_jeanne() = \t%lu nanoseconds\n", diff2 / nb_iteration);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
    for (i=0; i<nb_iteration; i++)
    {
        mult_jeanne_neon(res, c, d, f, STACK_1, gen_B, gen_Mtiprime16, gen_Mti16, rec);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
    diff3 = BILLION * (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec);

    printf("mean time using mult_jeanne_neon() = \t%lu nanoseconds\n", diff3 / nb_iteration);

    if(!amns_param) {
        printf("Erreur lors de l'exécution sage\n");
    }

    return 0;
}