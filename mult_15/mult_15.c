#include "mult_15.h"

//--------------------------------------------------------------------------
//--------------------------- GENERAL --------------------------------------
//--------------------------------------------------------------------------

static inline amns_block mul_lambda(amns_block x, amns_block lambda) {
    
    switch(lambda) {
        case 0:     return 0;

        case 1:     return x;
        case -1:    return -x;

        case 2:     return x << 1;
        case -2:    return -(x << 1);
        
        case 3:     return (x << 1) + x;
        case -3:    return -((x << 1) + x);

        case 4:     return x << 2;
        case -4:    return -(x << 2);

        case 5:     return (x << 2) + x;
        case -5:    return -((x << 2) + x);

        case 6:     return (x << 2) + (x << 1);
        case -6:    return -((x << 2) + (x << 1));

        case 7:     return (x << 3) - x;
        case -7:    return -((x << 3) - x);

        case 8:     return x << 3;
        case -8:    return -(x << 3);

        case 9:     return (x << 3) + x;
        case -9:    return -((x << 3) + x);

        case 10:     return x << 4;
        case -10:    return -(x << 4);

        default:
            return x * lambda;
    }
}

void get_gen_from_vec(amns_block gen[29], amns_block B[15], amns_block lambda) {
    /*
    Construit un vecteur représentant la partie génératrice de la matrice de Toeplitz associée au vecteur B (sans construire la matrice).
    La matrice associée aurait été de la forme :
    [
        B,
        x*B mod E,
        (x^2)*B mod E,
        ...
        (x^(n-1))*B mod E
    ]

    avec E(X) = X^15 - lambda.
    Ici on ne récupère que les éléments de la première ligne et de la première colonne (du coin inférieur gauche au coin supérieur droit).
    */

    for (int i = 0; i < 14; i++) {
        gen[i] = mul_lambda(B[1 + i], lambda);
    }

    for (int i = 0; i < 15; i++) {
        gen[14 + i] = B[i];
    }
}

void get_gen_15(amns_block gen_Mti[29], amns_elt * Mti) {
    /*
    Construit un vecteur correspondant à la partie "génératrice" de la matrice de Toeplitz Mti.
    On commence en haut à droite et finit en bas à gauche (les matrices de Toeplitz n'ont pas été construites à ma façon mais de façon transposée).
    */
    
    for (int i = 0; i < 15; i++) {
        gen_Mti[i] = Mti[0][14 - i];
    }
    
    for (int i = 1; i < 15; i++) {
        gen_Mti[14 + i] = Mti[i][0];
    }
}

static inline void mult_5_M(amns_llong * restrict result, const amns_llong * restrict V, const amns_block * restrict M) {
    /*
    Multiplie un vecteur V de longueur 5 avec un vecteur M de longueur 10 (c'est un générateur de matrice de Toeplitz).
    */

    for (int i = 0; i < 5; i++) {
        result[i] = 0;
        for (int j = 0; j < 5; j++) {
            result[i] += V[j] * M[4 - j + i];
        }
    }
}

//--------------------------------------------------------------------------
//--------------------------- SANS NEON ------------------------------------
//--------------------------------------------------------------------------

static inline void mult_5(amns_llong * restrict result, const amns_block * restrict V, const amns_block * restrict M) {
    /*
    Multiplie un vecteur V de longueur 5 avec un vecteur M de longueur 10 (c'est un générateur de matrice de Toeplitz).
    */
    
    for (int i = 0; i < 5; i++) {
        result[i] = 0;
        for (int j = 0; j < 5; j++) {
            result[i] += (amns_llong)V[j] * M[4 - j + i];
        }
    }
}

static void mult_toeplitz_15_gen(amns_llong * restrict result, amns_block * restrict V, amns_block * restrict gen_B) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 15*15.

    Note : A*B mod E plus perfomant en Toeplitz (50/60 ns de moins que schoolbook).
    */

    amns_block *gen = gen_B;

    // for(int i = 0; i < 15; i++)
    // {
    //     result[i] = 0;
    //     for(int j = 0; j < 15; j++)
    //         result[i] += (int64_t)V[j]*gen[14-j+i];
    // }

    amns_llong p0[5], p1[5], p2[5], p3[5], p4[5], p5[5];
    amns_block v0m2[5], v1m2[5], v0m1[5];
    const amns_block *v1 = V + 5;
    const amns_block *v2 = V + 10;
    amns_block m012[9], m123[9], m234[9];

    for (int i = 0; i < 5; i++) {
        v0m2[i] = V[i] - v2[i];
        v1m2[i] = v1[i] - v2[i];
        v0m1[i] = V[i] - v1[i];
    }

    for (int i = 0; i < 9; i++) {
        amns_block g0 = gen[i];
        amns_block g1 = gen[5 + i];
        amns_block g2 = gen[10 + i];
        amns_block g3 = gen[15 + i];
        amns_block g4 = gen[20 + i];

    	amns_block g5_10 = g1 + g2;
    	amns_block g10_15 = g2 + g3;

    	m012[i] = g0 + g5_10;
    	m123[i] = g5_10 + g3;
    	m234[i] = g10_15 + g4;
     }

    mult_5(p0, V + 10, m012);
    mult_5(p1, V + 5, m123);
    mult_5(p2, V, m234);

    mult_5(p3, v0m2, gen + 10);
    mult_5(p4, v1m2, gen + 5);
    mult_5(p5, v0m1, gen + 15);

    for (int i = 0; i < 5; i++) {
        result[i] = p0[i] + p3[i] + p4[i];
        result[i + 5] = p1[i] + p5[i] - p4[i];
        result[i + 10] = p2[i] - p3[i] - p5[i];
    }
}

static inline void mult_5_Mp(amns_block * restrict result, const amns_block * restrict V, const amns_block * restrict M) {
    /*
    Multiplie un vecteur V de longueur 5 avec un vecteur M de longueur 10 (c'est un générateur de matrice de Toeplitz).
    */

    for (int i = 0; i < 5; i++) {
        result[i] = 0;
        for (int j = 0; j < 5; j++) {
            result[i] += V[j] * M[4 - j + i];
        }
    }
}

static void mult_toeplitz_15_Mp(amns_block * restrict result, amns_llong * restrict V, amns_block gen_Mtiprime15[29]) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 15*15.

    Note : Plus performant en Toeplitz (20/30 ns de moins que schoolbook).
    */

    amns_block *gen_Mti = gen_Mtiprime15;

    // for(int i = 0; i < 15; i++)
    // {
    //     result[i] = 0;
    //     for(int j = 0; j < 15; j++)
    //         result[i] += V[j]*gen_Mti[14-j+i];
    // }

    amns_block v[15];

    for(int i = 0; i < 15; i++) {
        v[i] = V[i];
    }

    amns_block p0[5], p1[5], p2[5], p3[5], p4[5], p5[5];
    amns_block v0m2[5], v1m2[5], v0m1[5];
    const amns_block *v1 = v + 5;
    const amns_block *v2 = v + 10;
    amns_block m03, m034[9], m013[9], m012[9];

    for (int i = 0; i < 5; i++) {
        v0m2[i] = v[i] - v2[i];
        v1m2[i] = v1[i] - v2[i];
        v0m1[i] = v[i] - v1[i];
    }

    for (int i = 0; i < 9; i++) {
        m03 = gen_Mti[i + 10] + gen_Mti[i + 5];
        m034[i] = m03 + gen_Mti[i];
        m013[i] = m03 + gen_Mti[i + 15];
        m012[i] = gen_Mti[10 + i] + gen_Mti[15 + i] + gen_Mti[20 + i];
    }

    mult_5_Mp(p0, v + 10, m034);
    mult_5_Mp(p1, v + 5, m013);
    mult_5_Mp(p2, v, m012);

    mult_5_Mp(p3, v1m2, gen_Mti + 5);
    mult_5_Mp(p4, v0m2, gen_Mti + 10);
    mult_5_Mp(p5, v0m1, gen_Mti + 15);

    for (int i = 0; i < 5; i++) {
        result[i] = p0[i] + p3[i] + p4[i];
        result[i + 5] = p1[i] - p3[i] + p5[i];
        result[i + 10] = p2[i] - p4[i] - p5[i];
    }
}

static void mult_toeplitz_15_M(amns_llong * restrict result, amns_block * restrict V, amns_block gen_Mti15[29]) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 15*15.

    Note : Plus performant en schoolbook (80/90 ns de moins que Toeplitz).
    */

    amns_block *gen_Mti = gen_Mti15;

    for(int i = 0; i < 15; i++)
    {
        result[i] = 0;
        for(int j = 0; j < 15; j++)
            result[i] += (int64_t)V[j]*gen_Mti[14-j+i];
    }

    // amns_llong p0[5], p1[5], p2[5], p3[5], p4[5], p5[5];
    // amns_llong v0m2[5], v1m2[5], v0m1[5];
    // const amns_block *v1 = V + 5;
    // const amns_block *v2 = V + 10;
    // amns_block m012[9], m123[9], m234[9];

    // for (int i = 0; i < 5; i++) {
    //     v0m2[i] = (int64_t)V[i] - v2[i];
    //     v1m2[i] = (int64_t)v1[i] - v2[i];
    //     v0m1[i] = (int64_t)V[i] - v1[i];
    // }

    // for (int i = 0; i < 9; i++) {
    //     amns_block g0 = gen_Mti[i];
    //     amns_block g1 = gen_Mti[5 + i];
    //     amns_block g2 = gen_Mti[10 + i];
    //     amns_block g3 = gen_Mti[15 + i];
    //     amns_block g4 = gen_Mti[20 + i];

    // 	amns_block g5_10 = g1 + g2;
    // 	amns_block g10_15 = g2 + g3;

    // 	m012[i] = g0 + g5_10;
    // 	m123[i] = g5_10 + g3;
    // 	m234[i] = g10_15 + g4;
    // }

    // mult_5(p0, V + 10, m012);
    // mult_5(p1, V + 5, m123);
    // mult_5(p2, V, m234);

    // mult_5_M(p3, v0m2, gen_Mti + 10);
    // mult_5_M(p4, v1m2, gen_Mti + 5);
    // mult_5_M(p5, v0m1, gen_Mti + 15);

    // for (int i = 0; i < 5; i++) {
    //     result[i] = p0[i] + p3[i] + p4[i];
    //     result[i + 5] = p1[i] + p5[i] - p4[i];
    //     result[i + 10] = p2[i] - p3[i] - p5[i];
    // }
}

static void amns_internal_red_15_jeanne(amns_elt_ptr S, amns_llong * V, amns_srcptr AMNS, amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]) {
    /* 
    Réduit les coefficients de V pour qu'il reste dans l'AMNS.
    */
    
    int8_t i;

    amns_block Q[15];
    mult_toeplitz_15_Mp(Q, V, gen_Mtiprime15);

    amns_llong T[15];
    mult_toeplitz_15_M(T, Q, gen_Mti15);

    for (i=0; i<15; i++) {
        amns_llong tmp = V[i] + T[i];
        S[i] = (amns_block)(tmp >> AMNS_WORD_SIZE);
    }
}

static void amns_elt_mul_jeanne(amns_elt_ptr restrict c, amns_elt_ptr restrict a, amns_elt_ptr restrict b, amns_srcptr restrict AMNS, amns_block gen_B[29], amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]) {

    amns_llong tm[15] __attribute__((aligned(64)));
    mult_toeplitz_15_gen(tm, a, gen_B);

    amns_internal_red_15_jeanne(c, tm, AMNS, gen_Mtiprime15, gen_Mti15);
}

static void fp_elt_mul_jeanne(fp_elt_ptr dst, fp_elt_srcptr src1, fp_elt_srcptr src2, const fp_param param, uint8_t stack, amns_block gen_B[29], amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]) {

    amns_elt_mul_jeanne(*dst, *src1, *src2, param->AMNS, gen_B, gen_Mtiprime15, gen_Mti15);
}

void mult_jeanne(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[29], amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]) {

    fp_elt_mul_jeanne((fp_elt_ptr)(*c), (fp_elt_srcptr)(*a), (fp_elt_srcptr)(*b), f->param, stack, gen_B, gen_Mtiprime15, gen_Mti15);
}

//--------------------------------------------------------------------------
//--------------------------- AVEC NEON ------------------------------------
//--------------------------------------------------------------------------

static inline void mult_5_neon(amns_llong * restrict result, const amns_block * restrict V, const amns_block * restrict M) {
    /*
    Multiplie un vecteur V de longueur 5 avec un vecteur M de longueur 10 (c'est un générateur de matrice de Toeplitz).

    Note : Très très mauvais, plutôt utiliser mult_5_M.
    */

    int64x2_t acc;
    int32x2_t vV;
    int32x2_t vM;
    int64_t tmp0, tmp1;
    int32_t V_pad[8] = {V[0], V[1], V[2], V[3], V[4], 0, 0, 0};

    for (int i = 0; i < 5; i++) {

        int64x2_t acc0 = vdupq_n_s64(0);
        int64x2_t acc1 = vdupq_n_s64(0);

        int32x2_t vV0 = vld1_s32(&V_pad[0]);
        int32x2_t vM0 = vld1_s32(&M[4 + i]);
        acc0 = vmlal_s32(acc0, vV0, vM0);

        int32x2_t vV1 = vld1_s32(&V_pad[2]);
        int32x2_t vM1 = vld1_s32(&M[2 + i]);
        acc1 = vmlal_s32(acc1, vV1, vM1);

        int64_t tmp = vgetq_lane_s64(acc0, 0) + vgetq_lane_s64(acc0, 1) + vgetq_lane_s64(acc1, 0) + vgetq_lane_s64(acc1, 1);

        tmp += (int64_t)V[4]*M[i];

        tmp0 = vgetq_lane_s64(acc, 0);
        tmp1 = vgetq_lane_s64(acc, 1);

        result[i] = tmp;
    }
}

static void mult_toeplitz_15_gen_neon(amns_llong * restrict result, amns_block * restrict V, amns_block * restrict gen_B) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 15*15.

    Note : A*B mod E plus perfomant en Toeplitz (50/60 ns de moins que schoolbook).
    */

    amns_block *gen = gen_B;

    // for(int i = 0; i < 15; i++)
    // {
    //     result[i] = 0;
    //     for(int j = 0; j < 15; j++)
    //         result[i] += (int64_t)V[j]*gen[14-j+i];
    // }

    amns_llong p0[5], p1[5], p2[5], p3[5], p4[5], p5[5];
    amns_block v0m2[5], v1m2[5], v0m1[5];
    const amns_block *v1 = V + 5;
    const amns_block *v2 = V + 10;
    amns_block m012[9], m123[9], m234[9];

    for (int i = 0; i < 5; i++) {
        v0m2[i] = V[i] - v2[i];
        v1m2[i] = v1[i] - v2[i];
        v0m1[i] = V[i] - v1[i];
    }

    for (int i = 0; i < 9; i++) {
        amns_block g0 = gen[i];
        amns_block g1 = gen[5 + i];
        amns_block g2 = gen[10 + i];
        amns_block g3 = gen[15 + i];
        amns_block g4 = gen[20 + i];

    	amns_block g5_10 = g1 + g2;
    	amns_block g10_15 = g2 + g3;

    	m012[i] = g0 + g5_10;
    	m123[i] = g5_10 + g3;
    	m234[i] = g10_15 + g4;
     }

    // mult_5_neon(p0, V + 10, m012);
    // mult_5_neon(p1, V + 5, m123);
    // mult_5_neon(p2, V, m234);

    // mult_5_neon(p3, v0m2, gen + 10);
    // mult_5_neon(p4, v1m2, gen + 5);
    // mult_5_neon(p5, v0m1, gen + 15);

    mult_5(p0, V + 10, m012);
    mult_5(p1, V + 5, m123);
    mult_5(p2, V, m234);

    mult_5(p3, v0m2, gen + 10);
    mult_5(p4, v1m2, gen + 5);
    mult_5(p5, v0m1, gen + 15);

    for (int i = 0; i < 5; i++) {
        result[i] = p0[i] + p3[i] + p4[i];
        result[i + 5] = p1[i] + p5[i] - p4[i];
        result[i + 10] = p2[i] - p3[i] - p5[i];
    }
}

static void mult_5_Mp_neon(amns_block * restrict result, const amns_block * restrict V, const amns_block * restrict M) {
    /*
    Multiplie un vecteur V de longueur 5 avec un vecteur M de longueur 10 (c'est un générateur de matrice de Toeplitz).
    */

    // --- Bloc 0 : result[0..3] ---
    int32x4_t acc0 = vdupq_n_s32(0);

    for (int j = 0; j < 5; j++) {
        int32x4_t vV = vdupq_n_s32(V[j]);
        int32x4_t vM = vld1q_s32(&M[4 - j]);
        acc0 = vmlaq_s32(acc0, vV, vM);
    }

    vst1q_s32(&result[0], acc0);

    // --- Bloc 1 : result[4] ---
    int32_t sum = 0;

    for (int j = 0; j < 5; j++) {
        sum += V[j]*M[4 - j + 4];
    }

    result[4] = sum;
}

static void mult_toeplitz_15_Mp_neon(amns_block * restrict result, amns_llong * restrict V, amns_block gen_Mtiprime15[29]) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 15*15.

    Note : Plus performant en Toeplitz (20/30 ns de moins que schoolbook).
    */

    amns_block *gen_Mti = gen_Mtiprime15;

    // for(int i = 0; i < 15; i++)
    // {
    //     result[i] = 0;
    //     for(int j = 0; j < 15; j++)
    //         result[i] += V[j]*gen_Mti[14-j+i];
    // }

    amns_block v[15];

    for(int i = 0; i < 15; i++) {
        v[i] = V[i];
    }

    amns_block p0[5], p1[5], p2[5], p3[5], p4[5], p5[5];
    amns_block v0m2[5], v1m2[5], v0m1[5];
    const amns_block *v1 = v + 5;
    const amns_block *v2 = v + 10;
    amns_block m03, m034[9], m013[9], m012[9];

    for (int i = 0; i < 5; i++) {
        v0m2[i] = v[i] - v2[i];
        v1m2[i] = v1[i] - v2[i];
        v0m1[i] = v[i] - v1[i];
    }

    for (int i = 0; i < 9; i++) {
        m03 = gen_Mti[i + 10] + gen_Mti[i + 5];
        m034[i] = m03 + gen_Mti[i];
        m013[i] = m03 + gen_Mti[i + 15];
        m012[i] = gen_Mti[10 + i] + gen_Mti[15 + i] + gen_Mti[20 + i];
    }

    mult_5_Mp_neon(p0, v + 10, m034);
    mult_5_Mp_neon(p1, v + 5, m013);
    mult_5_Mp_neon(p2, v, m012);

    mult_5_Mp_neon(p3, v1m2, gen_Mti + 5);
    mult_5_Mp_neon(p4, v0m2, gen_Mti + 10);
    mult_5_Mp_neon(p5, v0m1, gen_Mti + 15);

    for (int i = 0; i < 5; i++) {
        result[i] = p0[i] + p3[i] + p4[i];
        result[i + 5] = p1[i] - p3[i] + p5[i];
        result[i + 10] = p2[i] - p4[i] - p5[i];
    }
}

static void mult_toeplitz_15_M_neon(amns_llong * restrict result, amns_block * restrict V, amns_block gen_Mti15[29]) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 15*15.

    Note : Plus performant en schoolbook (80/90 ns de moins que Toeplitz).
    */

    amns_block *gen_Mti = gen_Mti15;

    int i;
    for(i = 0; i <= 11; i += 4) {
        int64x2_t acc0 = vdupq_n_s64(0);
        int64x2_t acc1 = vdupq_n_s64(0);

        for(int j = 0; j < 15; j++) {
            int32x4_t vM = vld1q_s32(&gen_Mti[14 + i - j]);

            int32x2_t vM_low = vget_low_s32(vM);
            int32x2_t vM_high = vget_high_s32(vM);

            int32_t v = V[j];

            int32x2_t vV = vdup_n_s32(v);

            acc0 = vmlal_s32(acc0, vV, vM_low);
            acc1 = vmlal_s32(acc1, vV, vM_high);
        }

        vst1q_s64(&result[i], acc0);
        vst1q_s64(&result[i + 2], acc1);
    }

    for(; i < 15; i++) {
        int64_t acc = 0;
        for(int j = 0; j < 15; j++) {
            acc += (int64_t) V[j]*gen_Mti[14 - j + i];
        }
        result[i] = acc;
    }

    // amns_llong p0[5], p1[5], p2[5], p3[5], p4[5], p5[5];
    // amns_llong v0m2[5], v1m2[5], v0m1[5];
    // const amns_block *v1 = V + 5;
    // const amns_block *v2 = V + 10;
    // amns_block m012[9], m123[9], m234[9];

    // for (int i = 0; i < 5; i++) {
    //     v0m2[i] = (int64_t)V[i] - v2[i];
    //     v1m2[i] = (int64_t)v1[i] - v2[i];
    //     v0m1[i] = (int64_t)V[i] - v1[i];
    // }

    // for (int i = 0; i < 9; i++) {
    //     amns_block g0 = gen_Mti[i];
    //     amns_block g1 = gen_Mti[5 + i];
    //     amns_block g2 = gen_Mti[10 + i];
    //     amns_block g3 = gen_Mti[15 + i];
    //     amns_block g4 = gen_Mti[20 + i];

    // 	amns_block g5_10 = g1 + g2;
    // 	amns_block g10_15 = g2 + g3;

    // 	m012[i] = g0 + g5_10;
    // 	m123[i] = g5_10 + g3;
    // 	m234[i] = g10_15 + g4;
    // }

    // mult_5(p0, V + 10, m012);
    // mult_5(p1, V + 5, m123);
    // mult_5(p2, V, m234);

    // mult_5_M(p3, v0m2, gen_Mti + 10);
    // mult_5_M(p4, v1m2, gen_Mti + 5);
    // mult_5_M(p5, v0m1, gen_Mti + 15);

    // for (int i = 0; i < 5; i++) {
    //     result[i] = p0[i] + p3[i] + p4[i];
    //     result[i + 5] = p1[i] + p5[i] - p4[i];
    //     result[i + 10] = p2[i] - p3[i] - p5[i];
    // }
}

static void amns_internal_red_15_jeanne_neon(amns_elt_ptr S, amns_llong * V, amns_srcptr AMNS, amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]) {
    /* 
    Réduit les coefficients de V pour qu'il reste dans l'AMNS.
    */
    int8_t i;

    amns_block Q[15];
    mult_toeplitz_15_Mp_neon(Q, V, gen_Mtiprime15);

    amns_llong T[15];
    mult_toeplitz_15_M_neon(T, Q, gen_Mti15);

    for (i=0; i<15; i++) {
        amns_llong tmp = V[i] + T[i];
        S[i] = (amns_block)(tmp >> AMNS_WORD_SIZE);
    }
}

static void amns_elt_mul_jeanne_neon(amns_elt_ptr restrict c, amns_elt_ptr restrict a, amns_elt_ptr restrict b, amns_srcptr restrict AMNS, amns_block gen_B[29], amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]) {

    amns_llong tm[15] __attribute__((aligned(64)));
    mult_toeplitz_15_gen_neon(tm, a, gen_B);

    amns_internal_red_15_jeanne_neon(c, tm, AMNS, gen_Mtiprime15, gen_Mti15);
}

static void fp_elt_mul_jeanne_neon(fp_elt_ptr dst, fp_elt_srcptr src1, fp_elt_srcptr src2, const fp_param param, uint8_t stack, amns_block gen_B[29], amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]) {

    amns_elt_mul_jeanne_neon(*dst, *src1, *src2, param->AMNS, gen_B, gen_Mtiprime15, gen_Mti15);
}

void mult_jeanne_neon(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[29], amns_block gen_Mtiprime15[29], amns_block gen_Mti15[29]) {

    fp_elt_mul_jeanne_neon((fp_elt_ptr)(*c), (fp_elt_srcptr)(*a), (fp_elt_srcptr)(*b), f->param, stack, gen_B, gen_Mtiprime15, gen_Mti15);
}