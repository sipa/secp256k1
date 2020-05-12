/**********************************************************************
 * Copyright (c) 2018-2020 Andrew Poelstra, Jonas Nick                *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_MODULE_SCHNORRSIG_MAIN_
#define _SECP256K1_MODULE_SCHNORRSIG_MAIN_

#include "include/secp256k1.h"
#include "include/secp256k1_schnorrsig.h"
#include "hash.h"

/* Initializes SHA256 with fixed midstate. This midstate was computed by applying
 * SHA256 to SHA256("BIP340/nonce")||SHA256("BIP340/nonce"). */
static void secp256k1_nonce_function_bip340_sha256_tagged(secp256k1_sha256 *sha) {
    secp256k1_sha256_initialize(sha);
    sha->s[0] = 0xa96e75cbul;
    sha->s[1] = 0x74f9f0acul;
    sha->s[2] = 0xc49e3c98ul;
    sha->s[3] = 0x202f99baul;
    sha->s[4] = 0x8946a616ul;
    sha->s[5] = 0x4accf415ul;
    sha->s[6] = 0x86e335c3ul;
    sha->s[7] = 0x48d0a072ul;

    sha->bytes = 64;
}

/* Initializes SHA256 with fixed midstate. This midstate was computed by applying
 * SHA256 to SHA256("BIP340/aux")||SHA256("BIP340/aux"). */
static void secp256k1_nonce_function_bip340_sha256_tagged_aux(secp256k1_sha256 *sha) {
    secp256k1_sha256_initialize(sha);
    sha->s[0] = 0x5d74a872ul;
    sha->s[1] = 0xd57064d4ul;
    sha->s[2] = 0x89495becul;
    sha->s[3] = 0x910f46f5ul;
    sha->s[4] = 0xcbc6fd3eul;
    sha->s[5] = 0xaf05d9d0ul;
    sha->s[6] = 0xcb781ce6ul;
    sha->s[7] = 0x062930acul;

    sha->bytes = 64;
}

/* algo16 argument for nonce_function_bip340 to derive the nonce exactly as stated in BIP-340
 * by using the correct tagged hash function. */
static const unsigned char bip340_algo16[16] = "BIP340/nonce\0\0\0\0";

static int nonce_function_bip340(unsigned char *nonce32, const unsigned char *msg32, const unsigned char *key32, const unsigned char *xonly_pk32, const unsigned char *algo16, void *data) {
    secp256k1_sha256 sha;
    unsigned char masked_key[32];
    int i;

    if (algo16 == NULL) {
        return 0;
    }

    if (data != NULL) {
        secp256k1_nonce_function_bip340_sha256_tagged_aux(&sha);
        secp256k1_sha256_write(&sha, data, 32);
        secp256k1_sha256_finalize(&sha, masked_key);
        for (i = 0; i < 32; i++) {
            masked_key[i] ^= key32[i];
        }
    }

    /* Tag the hash with algo16 which is important to avoid nonce reuse across
     * algorithms. If this nonce function is used in BIP-340 signing as defined
     * in the spec, an optimized tagging implementation is used. */
    if (memcmp(algo16, bip340_algo16, 16) == 0) {
        secp256k1_nonce_function_bip340_sha256_tagged(&sha);
    } else {
        int algo16_len = 16;
        /* Remove terminating null bytes */
        while (algo16_len > 0 && !algo16[algo16_len - 1]) {
            algo16_len--;
        }
        secp256k1_sha256_initialize_tagged(&sha, algo16, algo16_len);
    }

    /* Hash (masked-)key||pk||msg using the tagged hash as per the spec */
    if (data != NULL) {
        secp256k1_sha256_write(&sha, masked_key, 32);
    } else {
        secp256k1_sha256_write(&sha, key32, 32);
    }
    secp256k1_sha256_write(&sha, xonly_pk32, 32);
    secp256k1_sha256_write(&sha, msg32, 32);
    secp256k1_sha256_finalize(&sha, nonce32);
    return 1;
}

const secp256k1_nonce_function_hardened secp256k1_nonce_function_bip340 = nonce_function_bip340;

#endif
