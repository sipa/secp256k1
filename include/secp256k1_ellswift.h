#ifndef SECP256K1_ELLSWIFT_H
#define SECP256K1_ELLSWIFT_H

#include "secp256k1.h"
#include "secp256k1_extrakeys.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This module provides an implementation of the ElligatorSwift encoding
 * for secp256k1 x-only public keys. Given a uniformly random x-only key, this
 * produces a 64-byte encoding that is indistinguishable from uniformly
 * random bytes.
 *
 * ElligatorSwift is described in https://eprint.iacr.org/2022/759 by
 * Chavez-Saab, Rodriguez-Henriquez, and Tibouchi.
 *
 * Let f be the function from pairs of field elements to point X coordinates,
 * defined as follows (all operations modulo p = 2^256 - 2^32 - 977)
 * f(u,t):
 * - Let C = 0xa2d2ba93507f1df233770c2a797962cc61f6d15da14ecd47d8d27ae1cd5f852,
 *   a square root of -3.
 * - If u=0, set u=1 instead.
 * - If t=0, set t=1 instead.
 * - If u^3 + t^2 + 7 = 0, multiply t by 2.
 * - Let p = u^3 + t^2 + 7
 * - Let m = u^3 - t^2 + 7
 * - Let v = (C * m / p - 1) * u / 2
 * - Let w = p / (C * t * u)
 * - Let x1 = v
 * - Let x2 = -u - v
 * - Let x3 = u + w^2
 * - Return the first of [x3,x2,x1] that is an X coordinate on the curve
 *   (at least one of them is, for any inputs u and t).
 *
 * Then an ElligatorSwift encoding of x consists of the 32-byte big-endian
 * encodings of field elements u and t concatenated, where f(u,t) = x.
 * The encoding algorithm is described in the paper, and effectively picks a
 * uniformly random pair (u,t) among those which encode x.
 */

/* Construct a 64-byte ElligatorSwift encoding of a given x-only pubkey.
 *
 *  Returns: 1 when pubkey is valid.
 *  Args:    ctx:        pointer to a context object
 *  Out:     ell64:      pointer to a 64-byte array to be filled
 *  In:      rnd32:      pointer to 32 bytes of entropy (must be unpredictable)
 *           pubkey:     a pointer to a secp256k1_xonly_pubkey containing an
 *                       initialized public key
 *
 * This function runs in variable time.
 */
SECP256K1_API int secp256k1_ellswift_encode(
    const secp256k1_context* ctx,
    unsigned char *ell64,
    const unsigned char *rnd32,
    const secp256k1_xonly_pubkey *pubkey
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4);

/** Decode a 64-bytes ElligatorSwift encoded x-only public key.
 *
 *  Returns: always 1
 *  Args:    ctx:        pointer to a context object
 *  Out:     pubkey:     pointer to a secp256k1_xonly_pubkey that will be filled
 *  In:      ell64:      pointer to a 64-byte array to decode
 *
 * This function runs in variable time.
 */
SECP256K1_API int secp256k1_ellswift_decode(
    const secp256k1_context* ctx,
    secp256k1_xonly_pubkey *pubkey,
    const unsigned char *ell64
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3);

#ifdef __cplusplus
}
#endif

#endif /* SECP256K1_ELLSWIFT_H */
