/***********************************************************************
 * Copyright (c) 2022 Pieter Wuille                                    *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 ***********************************************************************/

#ifndef SECP256K1_MODULE_ELLSWIFT_MAIN_H
#define SECP256K1_MODULE_ELLSWIFT_MAIN_H

#include "../../../include/secp256k1.h"
#include "../../../include/secp256k1_ellswift.h"
#include "../../hash.h"

/** c1 = (sqrt(-3)-1)/2 */
static const secp256k1_fe secp256k1_ellswift_c1 = SECP256K1_FE_CONST(0x851695d4, 0x9a83f8ef, 0x919bb861, 0x53cbcb16, 0x630fb68a, 0xed0a766a, 0x3ec693d6, 0x8e6afa40);
/** c2 = (-sqrt(-3)-1)/2 = -(c1+1) */
static const secp256k1_fe secp256k1_ellswift_c2 = SECP256K1_FE_CONST(0x7ae96a2b, 0x657c0710, 0x6e64479e, 0xac3434e9, 0x9cf04975, 0x12f58995, 0xc1396c28, 0x719501ee);
/** c3 = (-sqrt(-3)+1)/2 = -c1 = c2+1 */
static const secp256k1_fe secp256k1_ellswift_c3 = SECP256K1_FE_CONST(0x7ae96a2b, 0x657c0710, 0x6e64479e, 0xac3434e9, 0x9cf04975, 0x12f58995, 0xc1396c28, 0x719501ef);

/** Decode ElligatorSwift encoding (u, t) to a fraction xn/xd representing a curve X coordinate. */
static void secp256k1_ellswift_xswiftec_frac_var(secp256k1_fe* xn, secp256k1_fe* xd, const secp256k1_fe* u, const secp256k1_fe* t) {
    /* The implemented algorithm is the following (all operations in GF(p)):
     *
     * - c0 = sqrt(-3) = 0xa2d2ba93507f1df233770c2a797962cc61f6d15da14ecd47d8d27ae1cd5f852
     * - If u=0, set u=1.
     * - If t=0, set t=1.
     * - If u^3+7+t^2 = 0, set t=2*t.
     * - Let X=(u^3+7-t^2)/(2*t)
     * - Let Y=(X+t)/(c0*u)
     * - If x3=u+4*Y^2 is a valid x coordinate, return x3.
     * - If x2=(-X/Y-u)/2 is a valid x coordinare, return x2.
     * - Return x1=(X/Y-u)/2 (which is now guaranteed to be a valid x coordinate).
     *
     * Introducing s=t^2, g=u^3+7, and simplifying x1=-(x2+u) we get:
     *
     * - ...
     * - Let s=t^2
     * - Let g=u^3+7
     * - If g+s=0, set t=2*t, s=4*s
     * - Let X=(g-s)/(2*t)
     * - Let Y=(X+t)/(c0*u) = (g+s)/(2*c0*t*u)
     * - If x3=u+4*Y^2 is a valid x coordinate, return x3.
     * - If x2=(-X/Y-u)/2 is a valid x coordinate, return it.
     * - Return x1=-(x2+u).
     *
     * Now substitute Y^2 = -(g+s)^2/(12*s*u^2) and X/Y = c0*u*(g-s)/(g+s)
     *
     * - ...
     * - If g+s=0, set s=4*s
     * - If x3=u-(g+s)^2/(3*s*u^2) is a valid x coordinate, return it.
     * - If x2=(-c0*u*(g-s)/(g+s)-u)/2 is a valid x coordinate, return it.
     * - Return x1=(c0*u*(g-s)/(g+s)-u)/2.
     *
     * Simplifying x2 using 2 additional constants:
     *
     * - c1 = (c0-1)/2 = 0x851695d49a83f8ef919bb86153cbcb16630fb68aed0a766a3ec693d68e6afa40
     * - c2 = (-c0-1)/2 = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
     * - ...
     * - If x2=u*(c1*s+c2*g)/(g+s) is a valid x coordinate, return it.
     * - ...
     *
     * Writing x3 as a fraction:
     *
     * - ...
     * - If x3=(3*s*u^3-(g+s)^2)/(3*s*u^2)
     * - ...

     * Overall, we get:
     *
     * - c1 = 0x851695d49a83f8ef919bb86153cbcb16630fb68aed0a766a3ec693d68e6afa40
     * - c2 = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
     * - If u=0, set u=1.
     * - If t=0, set s=1, else set s=t^2
     * - Let g=u^3+7
     * - If g+s=0, set s=4*s
     * - If x3=(3*s*u^3-(g+s)^2)/(3*s*u^2) is a valid x coordinate, return it.
     * - If x2=u*(c1*s+c2*g)/(g+s) is a valid x coordinate, return it.
     * - Return x1=-(x2+u)
     */
    secp256k1_fe u1, s, g, p, d, n, l;
    u1 = *u;
    if (EXPECT(secp256k1_fe_normalizes_to_zero_var(&u1), 0)) u1 = secp256k1_fe_one;
    secp256k1_fe_sqr(&s, t);
    if (EXPECT(secp256k1_fe_normalizes_to_zero_var(t), 0)) s = secp256k1_fe_one;
    secp256k1_fe_sqr(&l, &u1); /* l = u^2 */
    secp256k1_fe_mul(&g, &l, &u1); /* g = u^3 */
    secp256k1_fe_add(&g, &secp256k1_fe_const_b); /* g = u^3 + 7 */
    p = g; /* p = g */
    secp256k1_fe_add(&p, &s); /* p = g+s */
    if (EXPECT(secp256k1_fe_normalizes_to_zero_var(&p), 0)) {
        secp256k1_fe_mul_int(&s, 4); /* s = 4*s */
        /* recompute p = g+s */
        p = g; /* p = g */
        secp256k1_fe_add(&p, &s); /* p = g+s */
    }
    secp256k1_fe_mul(&d, &s, &l); /* d = s*u^2 */
    secp256k1_fe_mul_int(&d, 3); /* d = 3*s*u^2 */
    secp256k1_fe_sqr(&l, &p); /* l = (g+s)^2 */
    secp256k1_fe_negate(&l, &l, 1); /* l = -(g+s)^2 */
    secp256k1_fe_mul(&n, &d, &u1); /* n = 3*s*u^3 */
    secp256k1_fe_add(&n, &l); /* n = 3*s*u^3-(g+s)^2 */
    if (secp256k1_ge_x_frac_on_curve_var(&n, &d)) {
        /* Return n/d = (3*s*u^3-(g+s)^2)/(3*s*u^2) */
        *xn = n;
        *xd = d;
        return;
    }
    *xd = p;
    secp256k1_fe_mul(&l, &secp256k1_ellswift_c1, &s); /* l = c1*s */
    secp256k1_fe_mul(&n, &secp256k1_ellswift_c2, &g); /* n = c2*g */
    secp256k1_fe_add(&n, &l); /* n = c1*s+c2*g */
    secp256k1_fe_mul(&n, &n, &u1); /* n = u*(c1*s+c2*g) */
    /* Possible optimization: in the invocation below, d^2 = (g+s)^2 is computed,
     * which we already have computed above. This could be deduplicated. */
    if (secp256k1_ge_x_frac_on_curve_var(&n, &p)) {
        /* Return n/p = u*(c1*s+c2*g)/(g+s) */
        *xn = n;
        return;
    }
    secp256k1_fe_mul(&l, &p, &u1); /* l = u*(g+s) */
    secp256k1_fe_add(&n, &l); /* n = u*(c1*s+c2*g)+u*g*s */
    secp256k1_fe_negate(xn, &n, 2); /* n = -u*(c1*s+c2*g)+u*g*s */
#ifdef VERIFY
    VERIFY_CHECK(secp256k1_ge_x_frac_on_curve_var(xn, &p));
#endif
    /* Return n/p = -(u*(c1*s+c2*g)/(g+s)+u) */
}

/** Decode ElligatorSwift encoding (u, t) to X coordinate. */
static void secp256k1_ellswift_xswiftec_var(secp256k1_fe* x, const secp256k1_fe* u, const secp256k1_fe* t) {
    secp256k1_fe xn, xd;
    secp256k1_ellswift_xswiftec_frac_var(&xn, &xd, u, t);
    secp256k1_fe_inv_var(&xd, &xd);
    secp256k1_fe_mul(x, &xn, &xd);
}

/** Decode ElligatorSwift encoding (u, t) to point P. */
static void secp256k1_ellswift_swiftec_var(secp256k1_ge* p, const secp256k1_fe* u, const secp256k1_fe* t) {
    secp256k1_fe x;
    secp256k1_ellswift_xswiftec_var(&x, u, t);
    secp256k1_ge_set_xo_var(p, &x, secp256k1_fe_is_odd(t));
}

/* Try to complete an ElligatorSwift encoding (u, t) for X coordinate x, given u and x.
 *
 * There may be up to 8 distinct t values such that (u, t) decodes back to x, but also
 * fewer, or none at all. Each such partial inverse can be accessed individually using a
 * distinct input argument c (in range 0-7), and some or all of these may return failure.
 * The following guarantees exist:
 * - Given (x, u), no two distinct c values give the same successful result t.
 * - Every successful result maps back to x through secp256k1_ellswift_xswiftec_var.
 * - Given (x, u), all t values that map back to x can be reached by combining the
 *   successful results from this function over all c values, with the exception of:
 *   - this function cannot be called with u=0
 *   - no result with t=0 will be returned
 *   - no result for which u^3 + t^2 + 7 = 0 will be returned.
 */
static int secp256k1_ellswift_xswiftec_inv_var(secp256k1_fe* t, const secp256k1_fe* x, const secp256k1_fe* u, int c) {
    /* The implemented algorithm is this (all arithmetic, except involving c, is mod p):
     *
     * - If (c & 2) = 0:
     *   - If (-x-u) is a valid X coordinate, fail.
     *   - If (c & 1) = 0, let v=x; otherwise, let v=-x-u.
     *   - Let s=-(u^3+7)/(u^2+u*v+v^2)
     * - If (c & 2) = 2:
     *   - Let s=x-u
     *   - If s=0, fail.
     *   - Let r=sqrt(-s*(4*(u^3+7)+3*u^2*s)); fail if it doesn't exist
     *   - If (c & 1) = 1:
     *     - If r=0, fail.
     *     - Let r=-r.
     *   - Let v=(r/s-u)/2.
     * - Let w=sqrt(s); fail if it doesn't exist.
     * - If (c & 4) = 4, let w=-w.
     * - Return w*(c1*u-v).
     */
    secp256k1_fe v = *x, um = *u; /* v = x; um = u */
    secp256k1_fe g, m, s, w;
    secp256k1_fe_normalize_weak(&v);
    secp256k1_fe_normalize_weak(&um);
    secp256k1_fe_sqr(&g, u); /* g = u^2 */
    secp256k1_fe_mul(&g, &g, u); /* g = u^3 */
    secp256k1_fe_add(&g, &secp256k1_fe_const_b); /* g = u^3+7 */
    if (!(c & 2)) {
        secp256k1_fe o;
        m = v; /* m = x */
        secp256k1_fe_add(&m, &um); /* m = x+u */
        secp256k1_fe_negate(&m, &m, 2); /* m = -(x+u) */
        if (secp256k1_ge_x_on_curve_var(&m)) return 0; /* test if -(x+u) on curve */
        if (c & 1) v = m; /* if c&1, v = -x-u */
        /* v = algorithm v */
        o = um; /* o = u */
        secp256k1_fe_add(&o, &v); /* o = u+v */
        secp256k1_fe_sqr(&o, &o); /* o = (u+v)^2 */
        secp256k1_fe_negate(&o, &o, 1); /* o = -(u+v)^2 */
        secp256k1_fe_mul(&s, &um, &v); /* s = u*v */
        secp256k1_fe_add(&s, &o); /* s = u*v-(u+v)^2 = -(u^2+u*v+v^2) */
        secp256k1_fe_mul(&m, &s, &g); /* m = -(u^3+7)*(u^2+u*v+v^2) [= algorithm s*(u^2+u*v+v^2)^2] */
        if (secp256k1_fe_jacobi_var(&m) < 0) return 0; /* algorithm s is square iff m is square */
        secp256k1_fe_inv_var(&s, &s); /* s = -1/(u^2+u*v+v^2) */
        secp256k1_fe_mul(&s, &s, &g); /* s = -(u^3+7)/(u^2+u*v+v^2) [= algorithm s] */
    } else {
        secp256k1_fe r2, r;
        secp256k1_fe_negate(&m, &um, 1); /* m = -u */
        s = m;
        secp256k1_fe_add(&s, &v); /* s = x-u [= algorithm s] */
        if (secp256k1_fe_normalizes_to_zero_var(&s)) return 0; /* test s=0 */
        if (secp256k1_fe_jacobi_var(&s) < 0) return 0; /* early squareness check of s */
        secp256k1_fe_normalize_weak(&g);
        secp256k1_fe_mul_int(&g, 4); /* g = 4*(u^3+7) */
        secp256k1_fe_sqr(&r2, &um); /* r2 = u^2 */
        secp256k1_fe_mul_int(&r2, 3); /* r2 = 3*u^2 */
        secp256k1_fe_mul(&r2, &r2, &s); /* r2 = 3*u^2*s */
        secp256k1_fe_add(&r2, &g); /* r2 = 4*(u^3+7)+3*u^2*s */
        secp256k1_fe_mul(&r2, &r2, &s); /* r2 = s*(4*(u^3+7)+3*u^2*s) */
        secp256k1_fe_negate(&r2, &r2, 1); /* r2 = -s*(4*(u^3+7)+3*u^2*s) */
        if (secp256k1_fe_jacobi_var(&r2) < 0) return 0;
        VERIFY_CHECK(secp256k1_fe_sqrt(&r, &r2)); /* r = sqrt(r2) [= algorithm r] */
        if (c & 1) {
            if (secp256k1_fe_normalizes_to_zero_var(&r)) return 0; /* test r=0 */
            secp256k1_fe_negate(&r, &r, 1); /* r=-r [= algorithm r] */
        }
        secp256k1_fe_inv_var(&v, &s); /* v = 1/s */
        secp256k1_fe_mul(&v, &v, &r); /* v = r/s */
        secp256k1_fe_add(&v, &m); /* v = r/s-u */
        secp256k1_fe_half(&v); /* v = (r/s-u)/2 [= algorithm v] */
    }
    VERIFY_CHECK(secp256k1_fe_sqrt(&w, &s)); /* w = sqrt(s) [= algorithm w] */
    if (!(c & 4)) secp256k1_fe_negate(&w, &w, 1); /* w = -w [= algorithm -w] */
    secp256k1_fe_mul(&um, &um, &secp256k1_ellswift_c3); /* um = c3*u = -c1*u */
    secp256k1_fe_add(&um, &v); /* um = v-c1*u */
    secp256k1_fe_mul(t, &w, &um); /* t = -w*(v-c1*u) = w*(c1-u) */
    return 1;
}

/** Find an ElligatorSwift encoding (u, t) for X coordinate x.
 *
 * hasher is a SHA256 object which a incrementing 4-byte counter is added to to
 * generate randomness for the rejection sampling in this function. Its size plus
 * 4 (for the counter) plus 9 (for the SHA256 padding) must be a multiple of 64
 * for efficiency reasons.
 */
static void secp256k1_ellswift_xelligatorswift_var(secp256k1_fe* u, secp256k1_fe* t, const secp256k1_fe* x, const secp256k1_sha256* hasher) {
    /* Pool of 3-bit branch values. */
    unsigned char branch_hash[32];
    /* Number of 3-bit values in branch_hash left. */
    int branches_left = 0;
    /* Field elements u and branch values are extracted from
     * SHA256(hasher || cnt) for consecutive values of cnt. cnt==0
     * is first used to populate a pool of 64 4-bit branch values. The 64 cnt
     * values that follow are used to generate field elements u. cnt==65 (and
     * multiples thereof) are used to repopulate the pool and start over, if
     * that were ever necessary. */
    uint32_t cnt = 0;
    VERIFY_CHECK((hasher->bytes + 4 + 9) % 64 == 0);
    while (1) {
        int branch;
        /* If the pool of branch values is empty, populate it. */
        if (branches_left == 0) {
            secp256k1_sha256 hash = *hasher;
            unsigned char buf4[4];
            buf4[0] = cnt;
            buf4[1] = cnt >> 8;
            buf4[2] = cnt >> 16;
            buf4[3] = cnt >> 24;
            ++cnt;
            secp256k1_sha256_write(&hash, buf4, 4);
            secp256k1_sha256_finalize(&hash, branch_hash);
            branches_left = 64;
        }
        /* Take a 3-bit branch value from the branch pool (top bit is discarded). */
        --branches_left;
        branch = (branch_hash[branches_left >> 1] >> ((branches_left & 1) << 2)) & 7;
        /* Compute a new u value by hashing. */
        {
            secp256k1_sha256 hash = *hasher;
            unsigned char buf4[4];
            unsigned char u32[32];
            buf4[0] = cnt;
            buf4[1] = cnt >> 8;
            buf4[2] = cnt >> 16;
            buf4[3] = cnt >> 24;
            ++cnt;
            secp256k1_sha256_write(&hash, buf4, 4);
            secp256k1_sha256_finalize(&hash, u32);
            if (!secp256k1_fe_set_b32(u, u32)) continue;
            if (secp256k1_fe_is_zero(u)) continue;
        }
        /* Find a remainder t, and return it if found. */
        if (secp256k1_ellswift_xswiftec_inv_var(t, x, u, branch)) {
            secp256k1_fe_normalize_var(t);
            break;
        }
    }
}

/** Find an ElligatorSwift encoding (u, t) for point P. */
static void secp256k1_ellswift_elligatorswift_var(secp256k1_fe* u, secp256k1_fe* t, const secp256k1_ge* p, const secp256k1_sha256* hasher) {
    secp256k1_ellswift_xelligatorswift_var(u, t, &p->x, hasher);
    if (secp256k1_fe_is_odd(t) != secp256k1_fe_is_odd(&p->y)) {
        secp256k1_fe_negate(t, t, 1);
        secp256k1_fe_normalize_var(t);
    }
}

int secp256k1_ellswift_encode(const secp256k1_context* ctx, unsigned char *ell64, const secp256k1_pubkey *pubkey, const unsigned char *rnd32) {
    secp256k1_ge p;
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(ell64 != NULL);
    ARG_CHECK(pubkey != NULL);
    ARG_CHECK(rnd32 != NULL);

    if (secp256k1_pubkey_load(ctx, &p, pubkey)) {
        static const unsigned char PREFIX[128 - 9 - 4 - 32 - 33] = "secp256k1_ellswift_encode";
        secp256k1_fe u, t;
        unsigned char p33[33];
        secp256k1_sha256 hash;

        /* Set up hasher state */
        secp256k1_sha256_initialize(&hash);
        secp256k1_sha256_write(&hash, PREFIX, sizeof(PREFIX));
        secp256k1_sha256_write(&hash, rnd32, 32);
        secp256k1_fe_get_b32(p33, &p.x);
        p33[32] = secp256k1_fe_is_odd(&p.y);
        secp256k1_sha256_write(&hash, p33, sizeof(p33));
        VERIFY_CHECK(hash.bytes == 128 - 9 - 4);

        /* Compute ElligatorSwift encoding and construct output. */
        secp256k1_ellswift_elligatorswift_var(&u, &t, &p, &hash);
        secp256k1_fe_get_b32(ell64, &u);
        secp256k1_fe_get_b32(ell64 + 32, &t);
        return 1;
    }
    /* Only returned in case the provided pubkey is invalid. */
    return 0;
}

int secp256k1_ellswift_create(const secp256k1_context* ctx, unsigned char *ell64, const unsigned char *seckey32, const unsigned char *rnd32) {
    secp256k1_ge p;
    secp256k1_fe u, t;
    secp256k1_sha256 hash;
    secp256k1_scalar seckey_scalar;
    static const unsigned char PREFIX[32] = "secp256k1_ellswift_create";
    static const unsigned char ZERO[32] = {0};
    int ret = 0;

    /* Sanity check inputs. */
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(ell64 != NULL);
    memset(ell64, 0, 64);
    ARG_CHECK(secp256k1_ecmult_gen_context_is_built(&ctx->ecmult_gen_ctx));
    ARG_CHECK(seckey32 != NULL);

    /* Compute (affine) public key */
    ret = secp256k1_ec_pubkey_create_helper(&ctx->ecmult_gen_ctx, &seckey_scalar, &p, seckey32);
    secp256k1_declassify(ctx, &p, sizeof(p)); /* not constant time in produced pubkey */
    secp256k1_fe_normalize_var(&p.x);
    secp256k1_fe_normalize_var(&p.y);

    /* Set up hasher state */
    secp256k1_sha256_initialize(&hash);
    secp256k1_sha256_write(&hash, PREFIX, sizeof(PREFIX));
    secp256k1_sha256_write(&hash, seckey32, 32);
    secp256k1_sha256_write(&hash, rnd32 ? rnd32 : ZERO, 32);
    secp256k1_sha256_write(&hash, ZERO, 32 - 9 - 4);
    secp256k1_declassify(ctx, &hash, sizeof(hash)); /* hasher gets to declassify private key */

    /* Compute ElligatorSwift encoding and construct output. */
    secp256k1_ellswift_elligatorswift_var(&u, &t, &p, &hash);
    secp256k1_fe_get_b32(ell64, &u);
    secp256k1_fe_get_b32(ell64 + 32, &t);

    secp256k1_memczero(ell64, 64, !ret);
    secp256k1_scalar_clear(&seckey_scalar);

    return ret;
}

int secp256k1_ellswift_decode(const secp256k1_context* ctx, secp256k1_pubkey *pubkey, const unsigned char *ell64) {
    secp256k1_fe u, t;
    secp256k1_ge p;
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(pubkey != NULL);
    ARG_CHECK(ell64 != NULL);

    secp256k1_fe_set_b32(&u, ell64);
    secp256k1_fe_normalize_var(&u);
    secp256k1_fe_set_b32(&t, ell64 + 32);
    secp256k1_fe_normalize_var(&t);
    secp256k1_ellswift_swiftec_var(&p, &u, &t);
    secp256k1_pubkey_save(pubkey, &p);
    return 1;
}

static int ellswift_xdh_hash_function_sha256(unsigned char *output, const unsigned char *x32, const unsigned char *ours64, const unsigned char *theirs64, void *data) {
    secp256k1_sha256 sha;

    (void)data;

    secp256k1_sha256_initialize(&sha);
    if (secp256k1_memcmp_var(ours64, theirs64, 64) <= 0) {
        secp256k1_sha256_write(&sha, ours64, 64);
        secp256k1_sha256_write(&sha, theirs64, 64);
    } else {
        secp256k1_sha256_write(&sha, theirs64, 64);
        secp256k1_sha256_write(&sha, ours64, 64);
    }
    secp256k1_sha256_write(&sha, x32, 32);
    secp256k1_sha256_finalize(&sha, output);

    return 1;
}

const secp256k1_ellswift_xdh_hash_function secp256k1_ellswift_xdh_hash_function_sha256 = ellswift_xdh_hash_function_sha256;
const secp256k1_ellswift_xdh_hash_function secp256k1_ellswift_xdh_hash_function_default = ellswift_xdh_hash_function_sha256;

int secp256k1_ellswift_xdh(const secp256k1_context* ctx, unsigned char *output, const unsigned char* theirs64, const unsigned char* ours64, const unsigned char* seckey32, secp256k1_ellswift_xdh_hash_function hashfp, void *data) {
    int ret = 0;
    int overflow;
    secp256k1_scalar s;
    secp256k1_fe xn, xd, px, u, t;
    unsigned char sx[32];

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(output != NULL);
    ARG_CHECK(theirs64 != NULL);
    ARG_CHECK(ours64 != NULL);
    ARG_CHECK(seckey32 != NULL);

    if (hashfp == NULL) {
        hashfp = secp256k1_ellswift_xdh_hash_function_default;
    }

    /* Load remote public key (as fraction). */
    secp256k1_fe_set_b32(&u, theirs64);
    secp256k1_fe_normalize_var(&u);
    secp256k1_fe_set_b32(&t, theirs64 + 32);
    secp256k1_fe_normalize_var(&t);
    secp256k1_ellswift_xswiftec_frac_var(&xn, &xd, &u, &t);

    /* Load private key (using one if invalid). */
    secp256k1_scalar_set_b32(&s, seckey32, &overflow);
    overflow = secp256k1_scalar_is_zero(&s);
    secp256k1_scalar_cmov(&s, &secp256k1_scalar_one, overflow);

    /* Compute shared X coordinate. */
    secp256k1_ecmult_const_xonly(&px, &xn, &xd, &s, 256, 1);
    secp256k1_fe_normalize(&px);
    secp256k1_fe_get_b32(sx, &px);

    /* Invoke hasher */
    ret = hashfp(output, sx, ours64, theirs64, data);

    memset(sx, 0, 32);
    secp256k1_fe_clear(&px);
    secp256k1_scalar_clear(&s);

    return !!ret & !overflow;
}

#endif
