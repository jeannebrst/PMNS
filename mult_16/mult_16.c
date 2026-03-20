#include "mult_16.h"

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

void get_gen_from_vec(amns_block gen[31], amns_block B[16], amns_block lambda) {
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

    avec E(X) = X^16 - lambda.
    Ici on ne récupère que les éléments de la première ligne et de la première colonne (du coin inférieur gauche au coin supérieur droit).
    */

    for (int i = 0; i < 15; i++) {
        gen[i] = mul_lambda(B[1 + i], lambda);
    }

    for (int i = 0; i < 16; i++) {
        gen[15 + i] = B[i];
    }
}

void get_gen_16(amns_block gen_Mti[31], amns_elt * Mti) {
    /*
    Construit un vecteur correspondant à la partie "génératrice" de la matrice de Toeplitz Mti.
    On commence en haut à droite et finit en bas à gauche (les matrices de Toeplitz n'ont pas été construites à ma façon mais de façon transposée).
    */
    
    for (int i = 0; i < 16; i++) {
        gen_Mti[i] = Mti[0][15 - i];
    }
    
    for (int i = 1; i < 16; i++) {
        gen_Mti[15 + i] = Mti[i][0];
    }
}

static inline void mult_4(amns_llong * restrict result, const amns_block * restrict V, const amns_block * restrict M) {
    /*
    Multiplie un vecteur V de longueur 4 avec un vecteur M de longueur 7 (c'est un générateur de matrice de Toeplitz) naïvement.
    */

    for (int i = 0; i < 4; i++) {
        result[i] = 0;
        for (int j = 0; j < 4; j++) {
            result[i] += (amns_llong)V[j] * M[3 - j + i];
        }
    }
}

static inline void mult_4_M(amns_llong * restrict result, const amns_llong * restrict V, const amns_block * restrict M) {
    /*
    Multiplie un vecteur V de longueur 4 avec un vecteur M de longueur 7 (c'est un générateur de matrice de Toeplitz) naïvement.
    */

    for (int i = 0; i < 4; i++) {
        result[i] = 0;
        for (int j = 0; j < 4; j++) {
            result[i] += V[j] * M[3 - j + i];
        }
    }
}

static void mult_toeplitz_8(amns_llong * restrict result, amns_block * restrict v, amns_block * restrict M) {
    /*
    Calcule un produit intermédiaire d'un vecteur V de longueur 8 avec un vecteur M de longueur 15 (c'est un générateur de matrice de Toeplitz).
    */

    amns_llong p0[4], p1[4], p2[4];
    amns_block m01[7], m12[7];
    amns_llong v0m1[4];
    const amns_block *v1 = v + 4;

    for (int i = 0; i < 4; i++) {
        v0m1[i] = (amns_llong)v[i] - v1[i];
    }

    for (int i = 0; i < 7; i++) {
        amns_block g0 = M[i];
        amns_block g1 = M[4 + i];
        amns_block g2 = M[8 + i];

        m01[i] = g0 + g1;
        m12[i] = g1 + g2;
    }

    mult_4(p0, v + 4, m01);
    mult_4(p1, v, m12);
    mult_4_M(p2, v0m1, M + 4);

    for (int i = 0; i < 4; i++) {
        result[i] = p0[i] + p2[i];
        result[i + 4] = p1[i] - p2[i];
    }
}

static void mult_toeplitz_8_M(amns_llong * restrict result, amns_llong * restrict v, amns_block * restrict M) {
    /*
    Calcule un produit intermédiaire d'un vecteur V de longueur 8 avec un vecteur M de longueur 15 (c'est un générateur de matrice de Toeplitz).
    */

    amns_llong p0[4], p1[4], p2[4];
    amns_llong v0p1[4];
    const amns_llong *v1 = v + 4;
    amns_block m0m1[7], m0m2[7];

    for (int i = 0; i < 4; i++) {
        v0p1[i] = (int64_t)v[i] + v1[i];
    }

    for (int i = 0; i < 7; i++) {
        amns_block g0 = M[i];
        amns_block g1 = M[4 + i];
        amns_block g2 = M[8 + i];

        m0m1[i] = g1 - g2;
        m0m2[i] = g1 - g0;
    }

    mult_4_M(p0, v0p1, M + 4);
    mult_4_M(p1, v, m0m1);
    mult_4_M(p2, v + 4, m0m2);

    for (int i = 0; i < 4; i++) {
        result[i] = p0[i] - p2[i];
        result[i + 4] = p0[i] - p1[i];
    }
}

//--------------------------------------------------------------------------
//--------------------------- SANS NEON ------------------------------------
//--------------------------------------------------------------------------

static void mult_toeplitz_16_gen(amns_llong * restrict result, amns_block * restrict V, amns_block * restrict gen_B) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 16*16.
    
    Note : A*B mod E plus perfomant en schoolbook (40/50 ns de moins que Toeplitz).
    */

    amns_block *gen = gen_B;

    for(int i = 0; i < 16; i++)
    {
        result[i] = 0;
        for(int j = 0; j < 16; j++)
            result[i] += (int64_t)V[j]*gen[15-j+i];
    }

    // int i;
    // for(i = 0; i < 16; i += 4) {
    //     int64x2_t acc0 = vdupq_n_s64(0);
    //     int64x2_t acc1 = vdupq_n_s64(0);

    //     for(int j = 0; j < 16; j++) {
    //         int32x4_t vM = vld1q_s32(&gen_Mti[15 + i - j]);

    //         int32x2_t vM_low = vget_low_s32(vM);
    //         int32x2_t vM_high = vget_high_s32(vM);

    //         int32_t v = V[j];

    //         int32x2_t vV = vdup_n_s32(v);

    //         acc0 = vmlal_s32(acc0, vV, vM_low);
    //         acc1 = vmlal_s32(acc1, vV, vM_high);
    //     }

    //     vst1q_s64(&result[i], acc0);
    //     vst1q_s64(&result[i + 2], acc1);
    // }

    // amns_llong p0[8], p1[8], p2[8];
    // amns_block v0m1[8];
    // const amns_block *v1 = V + 8;
    // amns_block m01[15], m12[15];

    // for (int i = 0; i < 8; i++) {
    //     v0m1[i] = V[i] - v1[i];
    // }

    // for (int i = 0; i < 15; i++) {
    //     amns_block g0 = gen[i];
    //     amns_block g1 = gen[8 + i];
    //     amns_block g2 = gen[16 + i];

    //     m01[i] = g0 + g1;
    //     m12[i] = g1 + g2;
    // }

    // mult_toeplitz_8(p0, V + 8, m01);
    // mult_toeplitz_8(p1, V, m12);
    // mult_toeplitz_8(p2, v0m1, gen + 8);

    // for (int i = 0; i < 8; i++) {
    //     result[i] = p0[i] + p2[i];
    //     result[i + 8] = p1[i] - p2[i];
    // }
}

static inline void mult_4_Mp(amns_block * restrict result, const amns_block * restrict V, const amns_block * restrict M) {
    /*
    Multiplie un vecteur V de longueur 4 avec un vecteur M de longueur 7 (c'est un générateur de matrice de Toeplitz) naïvement.
    */

    for (int i = 0; i < 4; i++) {
        result[i] = 0;
        for (int j = 0; j < 4; j++) {
            result[i] += (amns_llong)V[j] * M[3 - j + i];
        }
    }
}

static void mult_toeplitz_8_Mp(amns_block * restrict result, amns_block * restrict v, amns_block * restrict M) {
    /*
    Calcule un produit intermédiaire d'un vecteur V de longueur 8 avec un vecteur M de longueur 15 (c'est un générateur de matrice de Toeplitz).
    */

    amns_block p0[4], p1[4], p2[4];
    amns_block v0p1[4];
    const amns_block *v1 = v + 4;
    amns_block m0m1[7], m0m2[7];

    for (int i = 0; i < 4; i++) {
        v0p1[i] = v[i] + v1[i];
    }

    for (int i = 0; i < 7; i++) {
        amns_block g0 = M[i];
        amns_block g1 = M[4 + i];
        amns_block g2 = M[8 + i];

        m0m1[i] = g1 - g2;
        m0m2[i] = g1 - g0;
    }

    mult_4_Mp(p0, v0p1, M + 4);
    mult_4_Mp(p1, v, m0m1);
    mult_4_Mp(p2, v + 4, m0m2);

    for (int i = 0; i < 4; i++) {
        result[i] = p0[i] - p2[i];
        result[i + 4] = p0[i] - p1[i];
    }
}

static void mult_toeplitz_16_Mp(amns_block * restrict result, amns_llong * restrict V, amns_block gen_Mtiprime16[31]) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 16*16.
    
    Note : Toeplitz très légèrement plus performant que schoolbook.
    */

    amns_block *gen_Mti = gen_Mtiprime16;

    // for(int i = 0; i < 16; i++)
    // {
    //     result[i] = 0;
    //     for(int j = 0; j < 16; j++)
    //         result[i] += V[j]*gen_Mti[15-j+i];
    // }

    amns_block v[16];

    for(int i = 0; i < 16; i++) {
        v[i] = V[i];
    }

    amns_block p0[8], p1[8], p2[8];
    amns_block v0p1[8];
    const amns_block *v1 = v + 8;
    amns_block m0m1[15], m0m2[15];

    for (int i = 0; i < 8; i++) {
        v0p1[i] = v[i] + v1[i];
    }

    for (int i = 0; i < 15; i++) {
        amns_block g0 = gen_Mti[i];
        amns_block g1 = gen_Mti[8 + i];
        amns_block g2 = gen_Mti[16 + i];

        m0m1[i] = g1 - g2;
        m0m2[i] = g1 - g0;
    }

    mult_toeplitz_8_Mp(p0, v0p1, gen_Mti + 8);
    mult_toeplitz_8_Mp(p1, v, m0m1);
    mult_toeplitz_8_Mp(p2, v + 8, m0m2);

    for (int i = 0; i < 8; i++) {
        result[i] = p0[i] - p2[i];
        result[i + 8] = p0[i] - p1[i];
    }
}

static void mult_toeplitz_16_M(amns_llong * restrict result, amns_block * restrict V, amns_block gen_Mti16[31]) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 16*16.
    
    Note : Plus performant en schoolbook (100 ns de moins que Toeplitz).
    */

    amns_block *gen_Mti = gen_Mti16;

    for(int i = 0; i < 16; i++)
    {
        result[i] = 0;
        for(int j = 0; j < 16; j++)
            result[i] += (int64_t)V[j]*gen_Mti[15-j+i];
    }

    // amns_llong p0[8], p1[8], p2[8];
    // amns_llong v0p1[8];
    // const amns_block *v1 = V + 8;
    // amns_block m0m1[15], m0m2[15];

    // for (int i = 0; i < 8; i++) {
    //     v0p1[i] = (int64_t)V[i] + v1[i];
    // }

    // for (int i = 0; i < 15; i++) {
    //     amns_block g0 = gen_Mti[i];
    //     amns_block g1 = gen_Mti[8 + i];
    //     amns_block g2 = gen_Mti[16 + i];

    //     m0m1[i] = g1 - g2;
    //     m0m2[i] = g1 - g0;
    // }

    // mult_toeplitz_8_M(p0, v0p1, gen_Mti + 8);
    // mult_toeplitz_8(p1, V, m0m1);
    // mult_toeplitz_8(p2, V + 8, m0m2);

    // for (int i = 0; i < 8; i++) {
    //     result[i] = p0[i] - p2[i];
    //     result[i + 8] = p0[i] - p1[i];
    // }
}

static void amns_internal_red_16_jeanne(amns_elt_ptr S, amns_llong * V, amns_srcptr AMNS, amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]) {
    /* 
    Réduit les coefficients de V pour qu'il reste dans l'AMNS.
    */

    int8_t i;

    amns_block Q[16];
    mult_toeplitz_16_Mp(Q, V, gen_Mtiprime16);

    amns_llong T[16];
    mult_toeplitz_16_M(T, Q, gen_Mti16);

    for (i=0; i<16; i++) {
        amns_llong tmp = V[i] + T[i];
        S[i] = (amns_block)(tmp >> AMNS_WORD_SIZE);
    }
    
}

static void amns_elt_mul_jeanne(amns_elt_ptr restrict c, amns_elt_ptr restrict a, amns_elt_ptr restrict b, amns_srcptr restrict AMNS, amns_block gen_B[31], amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]) {

    amns_llong tm[16] __attribute__((aligned(64)));
    mult_toeplitz_16_gen(tm, a, gen_B);

    amns_internal_red_16_jeanne(c, tm, AMNS, gen_Mtiprime16, gen_Mti16);
}

static void fp_elt_mul_jeanne(fp_elt_ptr dst, fp_elt_srcptr src1, fp_elt_srcptr src2, const fp_param param, uint8_t stack, amns_block gen_B[31], amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]) {

    amns_elt_mul_jeanne(*dst, *src1, *src2, param->AMNS, gen_B, gen_Mtiprime16, gen_Mti16);
}

void mult_jeanne(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[31], amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]) {

    fp_elt_mul_jeanne((fp_elt_ptr)(*c), (fp_elt_srcptr)(*a), (fp_elt_srcptr)(*b), f->param, stack, gen_B, gen_Mtiprime16, gen_Mti16);
}

//--------------------------------------------------------------------------
//--------------------------- AVEC NEON ------------------------------------
//--------------------------------------------------------------------------

static void mult_toeplitz_16_gen_neon(amns_llong * restrict result, amns_block * restrict V, amns_block * restrict gen_B) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 16*16.
    
    Note : A*B mod E plus perfomant en schoolbook (40/50 ns de moins que Toeplitz).
    */

    amns_block *gen = gen_B;

    int i;
    for(i = 0; i < 16; i += 4) {
        int64x2_t acc0 = vdupq_n_s64(0);
        int64x2_t acc1 = vdupq_n_s64(0);

        for(int j = 0; j < 16; j++) {
            int32x4_t vM = vld1q_s32(&gen[15 + i - j]);

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

    // int i;
    // for(i = 0; i < 16; i += 4) {
    //     int64x2_t acc0 = vdupq_n_s64(0);
    //     int64x2_t acc1 = vdupq_n_s64(0);

    //     for(int j = 0; j < 16; j++) {
    //         int32x4_t vM = vld1q_s32(&gen_Mti[15 + i - j]);

    //         int32x2_t vM_low = vget_low_s32(vM);
    //         int32x2_t vM_high = vget_high_s32(vM);

    //         int32_t v = V[j];

    //         int32x2_t vV = vdup_n_s32(v);

    //         acc0 = vmlal_s32(acc0, vV, vM_low);
    //         acc1 = vmlal_s32(acc1, vV, vM_high);
    //     }

    //     vst1q_s64(&result[i], acc0);
    //     vst1q_s64(&result[i + 2], acc1);
    // }

    // amns_llong p0[8], p1[8], p2[8];
    // amns_block v0m1[8];
    // const amns_block *v1 = V + 8;
    // amns_block m01[15], m12[15];

    // for (int i = 0; i < 8; i++) {
    //     v0m1[i] = V[i] - v1[i];
    // }

    // for (int i = 0; i < 15; i++) {
    //     amns_block g0 = gen[i];
    //     amns_block g1 = gen[8 + i];
    //     amns_block g2 = gen[16 + i];

    //     m01[i] = g0 + g1;
    //     m12[i] = g1 + g2;
    // }

    // mult_toeplitz_8(p0, V + 8, m01);
    // mult_toeplitz_8(p1, V, m12);
    // mult_toeplitz_8(p2, v0m1, gen + 8);

    // for (int i = 0; i < 8; i++) {
    //     result[i] = p0[i] + p2[i];
    //     result[i + 8] = p1[i] - p2[i];
    // }
}

static inline void mult_4_Mp_neon(amns_block * restrict result, const amns_block * restrict V, const amns_block * restrict M) {
    /*
    Multiplie un vecteur V de longueur 4 avec un vecteur M de longueur 7 (c'est un générateur de matrice de Toeplitz).
    */

    // --- Bloc 0 : result[0..3] ---
    int32x4_t acc = vdupq_n_s32(0);

    acc = vmlaq_n_s32(acc, vld1q_s32(&M[3]), V[0]);
    acc = vmlaq_n_s32(acc, vld1q_s32(&M[2]), V[1]);
    acc = vmlaq_n_s32(acc, vld1q_s32(&M[1]), V[2]);
    acc = vmlaq_n_s32(acc, vld1q_s32(&M[0]), V[3]);

    vst1q_s32(result, acc);
}

static void mult_toeplitz_8_Mp_neon(amns_block * restrict result, amns_block * restrict v, amns_block * restrict M) {
    /*
    Calcule un produit intermédiaire d'un vecteur V de longueur 8 avec un vecteur M de longueur 15 (c'est un générateur de matrice de Toeplitz).
    */

    amns_block p0[4], p1[4], p2[4];
    amns_block v0p1[4];
    const amns_block *v1 = v + 4;
    amns_block m0m1[7], m0m2[7];

    for (int i = 0; i < 4; i++) {
        v0p1[i] = v[i] + v1[i];
    }

    for (int i = 0; i < 7; i++) {
        amns_block g0 = M[i];
        amns_block g1 = M[4 + i];
        amns_block g2 = M[8 + i];

        m0m1[i] = g1 - g2;
        m0m2[i] = g1 - g0;
    }

    mult_4_Mp_neon(p0, v0p1, M + 4);
    mult_4_Mp_neon(p1, v, m0m1);
    mult_4_Mp_neon(p2, v + 4, m0m2);

    for (int i = 0; i < 4; i++) {
        result[i] = p0[i] - p2[i];
        result[i + 4] = p0[i] - p1[i];
    }
}

static void mult_toeplitz_16_Mp_neon(amns_block * restrict result, amns_llong * restrict V, amns_block gen_Mtiprime16[31]) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 16*16.
    
    Note : Toeplitz très légèrement plus performant que schoolbook.
    */

    amns_block *gen_Mti = gen_Mtiprime16;

    // for(int i = 0; i < 16; i++)
    // {
    //     result[i] = 0;
    //     for(int j = 0; j < 16; j++)
    //         result[i] += V[j]*gen_Mti[15-j+i];
    // }

    amns_block v[16];

    for(int i = 0; i < 16; i++) {
        v[i] = V[i];
    }

    amns_block p0[8], p1[8], p2[8];
    amns_block v0p1[8];
    const amns_block *v1 = v + 8;
    amns_block m0m1[15], m0m2[15];

    for (int i = 0; i < 8; i++) {
        v0p1[i] = v[i] + v1[i];
    }

    for (int i = 0; i < 15; i++) {
        amns_block g0 = gen_Mti[i];
        amns_block g1 = gen_Mti[8 + i];
        amns_block g2 = gen_Mti[16 + i];

        m0m1[i] = g1 - g2;
        m0m2[i] = g1 - g0;
    }

    mult_toeplitz_8_Mp_neon(p0, v0p1, gen_Mti + 8);
    mult_toeplitz_8_Mp_neon(p1, v, m0m1);
    mult_toeplitz_8_Mp_neon(p2, v + 8, m0m2);

    for (int i = 0; i < 8; i++) {
        result[i] = p0[i] - p2[i];
        result[i + 8] = p0[i] - p1[i];
    }
}

static void mult_toeplitz_16_M_neon(amns_llong * restrict result, amns_block * restrict V, amns_block gen_Mti16[31]) {
    /*
    Calcule un produit vecteur-matrice dans le cas où la matrice est de Toeplitz et est 16*16.
    
    Note : Plus performant en schoolbook (100 ns de moins que Toeplitz).
    */

    amns_block *gen_Mti = gen_Mti16;

    int i;
    for(i = 0; i < 16; i += 4) {
        int64x2_t acc0 = vdupq_n_s64(0);
        int64x2_t acc1 = vdupq_n_s64(0);

        for(int j = 0; j < 16; j++) {
            int32x4_t vM = vld1q_s32(&gen_Mti[15 + i - j]);

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

    // amns_llong p0[8], p1[8], p2[8];
    // amns_llong v0p1[8];
    // const amns_block *v1 = V + 8;
    // amns_block m0m1[15], m0m2[15];

    // for (int i = 0; i < 8; i++) {
    //     v0p1[i] = (int64_t)V[i] + v1[i];
    // }

    // for (int i = 0; i < 15; i++) {
    //     amns_block g0 = gen_Mti[i];
    //     amns_block g1 = gen_Mti[8 + i];
    //     amns_block g2 = gen_Mti[16 + i];

    //     m0m1[i] = g1 - g2;
    //     m0m2[i] = g1 - g0;
    // }

    // mult_toeplitz_8_M(p0, v0p1, gen_Mti + 8);
    // mult_toeplitz_8(p1, V, m0m1);
    // mult_toeplitz_8(p2, V + 8, m0m2);

    // for (int i = 0; i < 8; i++) {
    //     result[i] = p0[i] - p2[i];
    //     result[i + 8] = p0[i] - p1[i];
    // }
}

static void amns_internal_red_16_jeanne_neon(amns_elt_ptr S, amns_llong * V, amns_srcptr AMNS, amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]) {
    /* 
    Réduit les coefficients de V pour qu'il reste dans l'AMNS.
    */

    int8_t i;

    amns_block Q[16];
    mult_toeplitz_16_Mp_neon(Q, V, gen_Mtiprime16);

    amns_llong T[16];
    mult_toeplitz_16_M_neon(T, Q, gen_Mti16);

    for (i=0; i<16; i++) {
        amns_llong tmp = V[i] + T[i];
        S[i] = (amns_block)(tmp >> AMNS_WORD_SIZE);
    }
    
}

static void amns_elt_mul_jeanne_neon(amns_elt_ptr restrict c, amns_elt_ptr restrict a, amns_elt_ptr restrict b, amns_srcptr restrict AMNS, amns_block gen_B[31], amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]) {

    amns_llong tm[16] __attribute__((aligned(64)));
    mult_toeplitz_16_gen_neon(tm, a, gen_B);

    amns_internal_red_16_jeanne_neon(c, tm, AMNS, gen_Mtiprime16, gen_Mti16);
}

static void fp_elt_mul_jeanne_neon(fp_elt_ptr dst, fp_elt_srcptr src1, fp_elt_srcptr src2, const fp_param param, uint8_t stack, amns_block gen_B[31], amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]) {

    amns_elt_mul_jeanne_neon(*dst, *src1, *src2, param->AMNS, gen_B, gen_Mtiprime16, gen_Mti16);
}

void mult_jeanne_neon(fe_ptr c, fe_srcptr a, fe_srcptr b, field_srcptr f, uint8_t stack, amns_block gen_B[31], amns_block gen_Mtiprime16[31], amns_block gen_Mti16[31]) {

    fp_elt_mul_jeanne_neon((fp_elt_ptr)(*c), (fp_elt_srcptr)(*a), (fp_elt_srcptr)(*b), f->param, stack, gen_B, gen_Mtiprime16, gen_Mti16);
}