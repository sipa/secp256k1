/***********************************************************************
 * Copyright (c) 2022 Pieter Wuille                                    *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 ***********************************************************************/

#ifndef SECP256K1_MODULE_ELLSWIFT_MAIN_H
#define SECP256K1_MODULE_ELLSWIFT_MAIN_H

#include "../../../include/secp256k1.h"
#include "../../../include/secp256k1_extrakeys.h"
#include "../../../include/secp256k1_ellswift.h"
#include "../../hash.h"

/* c1 = the square root of -3 ((-3)**((p+1)/4)). */
static const secp256k1_fe secp256k1_ellswift_c1 = SECP256K1_FE_CONST(0x0a2d2ba9, 0x3507f1df, 0x233770c2, 0xa797962c, 0xc61f6d15, 0xda14ecd4, 0x7d8d27ae, 0x1cd5f852);
/* c2 = -1/2 * (c1 - 1). */
static const secp256k1_fe secp256k1_ellswift_c2 = SECP256K1_FE_CONST(0x7ae96a2b, 0x657c0710, 0x6e64479e, 0xac3434e9, 0x9cf04975, 0x12f58995, 0xc1396c28, 0x719501ef);

/* Given field elements (u,t), compute a field element out which is a valid
 * X coordinate on the curve.
 * It implements the function f(u,t) defined in secp256k1_ellswift.h.
 */
static void secp256k1_ellswift_fe2_to_gex_var(secp256k1_fe* out, const secp256k1_fe* u, const secp256k1_fe* t) {
    secp256k1_fe v1 = *u, v2 = *t;
    secp256k1_fe v3, v4, v5, v6, v7, v8;
    secp256k1_fe_normalize_var(&v1);
    secp256k1_fe_normalize_var(&v2);
    if (secp256k1_fe_is_zero(&v1)) v1 = secp256k1_fe_one;
    if (secp256k1_fe_is_zero(&v2)) v2 = secp256k1_fe_one;
    secp256k1_fe_sqr(&v3, &v1);
    secp256k1_fe_mul(&v3, &v3, &v1);
    secp256k1_fe_add(&v3, &secp256k1_fe_const_b);
    secp256k1_fe_sqr(&v4, &v2);
    v5 = v3;
    secp256k1_fe_add(&v5, &v4);
    if (secp256k1_fe_normalizes_to_zero_var(&v5)) {
        secp256k1_fe_add(&v2, &v2);
        secp256k1_fe_sqr(&v4, &v2);
        v5 = v3;
        secp256k1_fe_add(&v5, &v4);
    }
    secp256k1_fe_mul(&v6, &v1, &secp256k1_ellswift_c1);
    secp256k1_fe_negate(&v4, &v4, 1);
    secp256k1_fe_add(&v4, &v3);
    secp256k1_fe_mul(&v4, &v4, &v6);
    secp256k1_fe_mul(&v2, &v2, &v6);
    secp256k1_fe_sqr(&v2, &v2);
    secp256k1_fe_sqr(&v8, &v5);
    secp256k1_fe_mul(&v3, &v1, &v2);
    secp256k1_fe_add(&v3, &v8);
    secp256k1_fe_sqr(&v6, &v2);
    secp256k1_fe_sqr(&v6, &v6);
    secp256k1_fe_mul_int(&v6, 7);
    secp256k1_fe_sqr(&v7, &v3);
    secp256k1_fe_mul(&v7, &v7, &v3);
    secp256k1_fe_mul(&v7, &v7, &v2);
    secp256k1_fe_add(&v7, &v6);
    if (secp256k1_fe_jacobi_var(&v7) >= 0) {
        secp256k1_fe_inv_var(out, &v2);
        secp256k1_fe_mul(out, out, &v3);
        return;
    }
    secp256k1_fe_mul(&v1, &v1, &v5);
    secp256k1_fe_add(&v1, &v4);
    secp256k1_fe_half(&v1);
    secp256k1_fe_negate(&v1, &v1, 3);
    secp256k1_fe_sqr(&v6, &v8);
    secp256k1_fe_mul_int(&v6, 7);
    secp256k1_fe_sqr(&v7, &v1);
    secp256k1_fe_mul(&v7, &v7, &v1);
    secp256k1_fe_mul(&v7, &v7, &v5);
    secp256k1_fe_add(&v7, &v6);
    secp256k1_fe_inv_var(&v5, &v5);
    if (secp256k1_fe_jacobi_var(&v7) >= 0) {
        secp256k1_fe_mul(out, &v5, &v1);
        return;
    }
    secp256k1_fe_add(&v1, &v4);
    secp256k1_fe_mul(out, &v5, &v1);
}

/* Given an X coordinate on the curve x, a non-zero field element u, and an
 * integer branch value i in [0, 8), compute a field element t (in out), such that
 * secp256k1_ellswift_fe2_to_gex_var(u, t) returns x, or fails. Combining all
 * non-failing outs for a given (x, u), over all values of i, results in the set
 * of all preimages of p under secp256k1_ellswift_fe2_to_gex_var. No
 * two (x, u, i) inputs map to the same out, if successful.
 */
static int secp256k1_ellswift_fegex_to_fe_var(secp256k1_fe* out, const secp256k1_fe* x, secp256k1_fe* u, int i) {
    secp256k1_fe xm = *x, um = *u;
    secp256k1_fe g, s, w2, w;
    secp256k1_fe_normalize_weak(&xm);
    secp256k1_fe_normalize_weak(&um);
    secp256k1_fe_sqr(&g, u);
    secp256k1_fe_mul(&g, &g, u);
    secp256k1_fe_add(&g, &secp256k1_fe_const_b);
    if ((i & 2) == 0) {
        secp256k1_fe o;
        s = xm;
        secp256k1_fe_add(&s, &um);
        secp256k1_fe_sqr(&o, &s);
        secp256k1_fe_mul(&o, &o, &s);
        secp256k1_fe_negate(&o, &o, 1);
        secp256k1_fe_add(&o, &secp256k1_fe_const_b);
        if (secp256k1_fe_jacobi_var(&o) >= 0) return 0;
        if (i & 1) {
            secp256k1_fe_add(&xm, &um);
            secp256k1_fe_negate(&xm, &xm, 2);
        }
        o = um;
        secp256k1_fe_add(&o, &xm);
        secp256k1_fe_sqr(&o, &o);
        secp256k1_fe_negate(&o, &o, 1);
        secp256k1_fe_mul(&w2, &um, &xm);
        secp256k1_fe_add(&w2, &o);
        secp256k1_fe_inv_var(&w2, &w2);
        secp256k1_fe_mul(&w2, &w2, &g);
    } else {
        secp256k1_fe r2, r;
        secp256k1_fe_negate(&w2, &um, 1);
        secp256k1_fe_add(&w2, &xm);
        if (secp256k1_fe_normalizes_to_zero_var(&w2)) return 0;
        secp256k1_fe_normalize_weak(&g);
        secp256k1_fe_mul_int(&g, 4);
        secp256k1_fe_sqr(&r2, &um);
        secp256k1_fe_mul_int(&r2, 3);
        secp256k1_fe_mul(&r2, &r2, &w2);
        secp256k1_fe_add(&r2, &g);
        secp256k1_fe_mul(&r2, &r2, &w2);
        secp256k1_fe_negate(&r2, &r2, 1);
        if (!secp256k1_fe_sqrt(&r, &r2)) return 0;
        if (i & 1) {
            if (secp256k1_fe_normalizes_to_zero_var(&r)) return 0;
            secp256k1_fe_negate(&r, &r, 1);
        }
        secp256k1_fe_inv_var(&xm, &w2);
        secp256k1_fe_mul(&xm, &xm, &r);
        secp256k1_fe_add(&xm, &um);
        secp256k1_fe_half(&xm);
        secp256k1_fe_negate(&xm, &xm, 2);
    }
    if (!secp256k1_fe_sqrt(&w, &w2)) return 0;
    if (i & 4) secp256k1_fe_negate(&w, &w, 1);
    secp256k1_fe_mul(&um, &um, &secp256k1_ellswift_c2);
    secp256k1_fe_add(&um, &xm);
    secp256k1_fe_mul(out, &w, &um);
    return 1;
}

int secp256k1_ellswift_encode(const secp256k1_context* ctx, unsigned char *ell64, const unsigned char *rnd32, const secp256k1_xonly_pubkey *pubkey) {
    secp256k1_ge p;
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(ell64 != NULL);
    ARG_CHECK(rnd32 != NULL);
    ARG_CHECK(pubkey != NULL);

    if (secp256k1_xonly_pubkey_load(ctx, &p, pubkey)) {
        uint32_t cnt = 0;
        /* Field elements and branch values are extracted from
         * SHA256("secp256k1_ellswift_encode\x00" + uint32{cnt} + rnd32 + X)
         * for consecutive values of cnt. cnt==0 is first used to populate a pool of
         * 64 4-bit branch values. The 64 cnt values that follow are used to
         * generate field elements u. cnt==65 (and multiples thereof) are used to
         * repopulate the pool and start over, if that were ever necessary. */
        unsigned char hashdata[26 + 4 + 32 + 32] = "secp256k1_ellswift_encode";
        /* Pool of 3-bit branch values. */
        unsigned char branch_hash[32];
        /* Number of 3-bit values in branch_hash left. */
        int branches_left = 0;
        /* Fill up hashdata, excluding i. */
        memcpy(hashdata + 26 + 4, rnd32, 32);
        secp256k1_fe_get_b32(hashdata + 26 + 4 + 32, &p.x);
        while (1) {
            int branch;
            secp256k1_fe u, t;
            /* If the pool of branch values is empty, populate it. */
            if (branches_left == 0) {
                secp256k1_sha256 hash;
                hashdata[26 + 0] = cnt;
                hashdata[26 + 1] = cnt >> 8;
                hashdata[26 + 2] = cnt >> 16;
                hashdata[26 + 3] = cnt >> 24;
                secp256k1_sha256_initialize(&hash);
                secp256k1_sha256_write(&hash, hashdata, sizeof(hashdata));
                secp256k1_sha256_finalize(&hash, branch_hash);
                ++cnt;
                branches_left = 64;
            }
            /* Take a 3-bit branch value from the branch pool (top bit is discarded). */
            --branches_left;
            branch = (branch_hash[(63 - branches_left) >> 1] >> (((63 - branches_left) & 1) << 2)) & 7;
            /* Compute a new u value by hashing (a potential first 32 bytes of the output). */
            {
                secp256k1_sha256 hash;
                hashdata[26 + 0] = cnt;
                hashdata[26 + 1] = cnt >> 8;
                hashdata[26 + 2] = cnt >> 16;
                hashdata[26 + 3] = cnt >> 24;
                secp256k1_sha256_initialize(&hash);
                secp256k1_sha256_write(&hash, hashdata, sizeof(hashdata));
                secp256k1_sha256_finalize(&hash, ell64);
                ++cnt;
            }
            if (!secp256k1_fe_set_b32(&u, ell64)) continue;
            if (secp256k1_fe_is_zero(&u)) continue;
            /* Compute the remainder t to encode in the last 32 bytes of the output. */
            if (secp256k1_ellswift_fegex_to_fe_var(&t, &p.x, &u, branch)) {
                secp256k1_fe_normalize_var(&t);
                secp256k1_fe_get_b32(ell64 + 32, &t);
                break;
            }
        }
        memset(hashdata, 0, sizeof(hashdata));
        return 1;
    }
    /* Only returned in case the provided pubkey is invalid. */
    return 0;
}

int secp256k1_ellswift_decode(const secp256k1_context* ctx, secp256k1_xonly_pubkey *pubkey, const unsigned char *ell64) {
    secp256k1_fe u, t;
    secp256k1_fe x;
    secp256k1_ge p;
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(pubkey != NULL);
    ARG_CHECK(ell64 != NULL);

    secp256k1_fe_set_b32(&u, ell64);
    secp256k1_fe_set_b32(&t, ell64 + 32);
    secp256k1_ellswift_fe2_to_gex_var(&x, &u, &t);
    secp256k1_ge_set_xo_var(&p, &x, 0);
    secp256k1_xonly_pubkey_save(pubkey, &p);
    return 1;
}

#endif
