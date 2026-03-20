from sage.all import *
import sys

load("amns_generator.sage")


#~ NOTE: Here, AMNS generation is done according to the implementation strategy mentioned at the end of section 5.1 of the article.

if __name__=="__main__":

    if len(sys.argv) != 3 :
        print("Usage Python : gen.py p n")
        sys.exit(1)

    word_size = 32

    p = Integer(int(sys.argv[1]))
    n = int(sys.argv[2])

    n_min = n
    n_max = n_min + 2

    abs_lamb_max = 1 << 3

    nth_root_max_duration_checks = 60  # seconds

    find_all_nthroot = True

    #~ amns generation
    flt_amns = generate_amns_candidates_with_n_min_max(word_size, p, n_min, n_max, abs_lamb_max, nth_root_max_duration_checks, find_all_nthroot)

    #~ Data structure for each AMNS generated: [delta, n, E, rho_log2, gamma, M, M', conv_P0, conv_P1]

    nb_amns = len(flt_amns)

    if nb_amns != 0:
        for cand in flt_amns:
            if cand[1] == n:
                print(cand)
                break