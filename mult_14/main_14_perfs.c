#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <arm_neon.h>

#include "mult_14.h"
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
    number_set_str(p, "4781668983906166242955001894344923773259119655253013193367", 10);

    // char *amns_param = get_amns_param("4781668983906166242955001894344923773259119655253013193367", "14");
    char *amns_param = "[11, 14, [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 19, 4490408646817400965549706395738165476182246170598781298672, [-435, 2568, 3941, 386, -4921, -6302, -2868, 4599, 3244, -5241, 914, -4171, 3602, -4066], [1588615239, 829992446, 1975879519, 1897731458, 3330844596, 2196031636, 2197856043, 3163498229, 2983258437, 4099481521, 537827624, 2363605116, 308238305, 2669670221], [922, -9879, -29366, -5526, -17677, 519, -2879, -11187, -8869, -5337, -5701, -10646, -4712, -2890], [-412, -9938, -9638, -13688, -3907, -2040, -1073, -18048, -4365, -7153, 2928, -9097, -8731, -2514]]";

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

    get_gen_14(gen_Mti14, ((fp_param)((field_srcptr)f)->param)->AMNS->Mti);
    get_gen_14(gen_Mtiprime14, ((fp_param)((field_srcptr)f)->param)->AMNS->Mtiprime);
    get_gen_from_vec(gen_B, *(fp_elt_srcptr)(*d), ((fp_param)((field_srcptr)f)->param)->AMNS->lambda);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
    for (i=0; i<nb_iteration; i++)
    {
        mult_jeanne(res, c, d, f, STACK_1, gen_B, gen_Mtiprime14, gen_Mti14);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
    diff2 = BILLION * (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec);

    printf("mean time using mult_jeanne() = \t%lu nanoseconds\n", diff2 / nb_iteration);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
    for (i=0; i<nb_iteration; i++)
    {
        mult_jeanne_neon(res, c, d, f, STACK_1, gen_B, gen_Mtiprime14, gen_Mti14);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
    diff3 = BILLION * (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec);

    printf("mean time using mult_jeanne_neon() = \t%lu nanoseconds\n", diff3 / nb_iteration);

    if(!amns_param) {
        printf("Erreur lors de l'exécution sage\n");
    }

    return 0;
}