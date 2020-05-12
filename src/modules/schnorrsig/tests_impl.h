/**********************************************************************
 * Copyright (c) 2018-2020 Andrew Poelstra, Jonas Nick                *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_MODULE_SCHNORRSIG_TESTS_

#define _SECP256K1_MODULE_SCHNORRSIG_TESTS_

#include "secp256k1_schnorrsig.h"

/* Checks that a bit flip in the n_flip-th argument (that has n_bytes many
 * bytes) changes the hash function
 */
void nonce_function_bip340_bitflip(unsigned char **args, size_t n_flip, size_t n_bytes) {
    unsigned char nonces[2][32];
    CHECK(nonce_function_bip340(nonces[0], args[0], args[1], args[2], args[3], args[4]) == 1);
    secp256k1_rand_flip(args[n_flip], n_bytes);
    CHECK(nonce_function_bip340(nonces[1], args[0], args[1], args[2], args[3], args[4]) == 1);
    CHECK(memcmp(nonces[0], nonces[1], 32) != 0);
}

void run_nonce_function_bip340_tests(void) {
    unsigned char tag[12] = "BIP340/nonce";
    unsigned char aux_tag[10] = "BIP340/aux";
    unsigned char algo16[16] = "BIP340/nonce\0\0\0\0";
    secp256k1_sha256 sha;
    secp256k1_sha256 sha_optimized;
    unsigned char nonce[32];
    unsigned char msg[32];
    unsigned char key[32];
    unsigned char pk[32];
    unsigned char aux_rand[32];
    unsigned char *args[5];
    int i;

    /* Check that hash initialized by
     * secp256k1_nonce_function_bip340_sha256_tagged has the expected
     * state. */
    secp256k1_sha256_initialize_tagged(&sha, tag, sizeof(tag));
    secp256k1_nonce_function_bip340_sha256_tagged(&sha_optimized);
    test_sha256_eq(&sha, &sha_optimized);

   /* Check that hash initialized by
    * secp256k1_nonce_function_bip340_sha256_tagged_aux has the expected
    * state. */
    secp256k1_sha256_initialize_tagged(&sha, aux_tag, sizeof(aux_tag));
    secp256k1_nonce_function_bip340_sha256_tagged_aux(&sha_optimized);
    test_sha256_eq(&sha, &sha_optimized);

    secp256k1_rand256(msg);
    secp256k1_rand256(key);
    secp256k1_rand256(pk);
    secp256k1_rand256(aux_rand);

    /* Check that a bitflip in an argument results in different nonces. */
    args[0] = msg;
    args[1] = key;
    args[2] = pk;
    args[3] = algo16;
    args[4] = aux_rand;
    for (i = 0; i < count; i++) {
        nonce_function_bip340_bitflip(args, 0, 32);
        nonce_function_bip340_bitflip(args, 1, 32);
        nonce_function_bip340_bitflip(args, 2, 32);
        /* Flip algo16 special case "BIP340/nonce" */
        nonce_function_bip340_bitflip(args, 3, 16);
        /* Flip algo16 again */
        nonce_function_bip340_bitflip(args, 3, 16);
        nonce_function_bip340_bitflip(args, 4, 32);
    }

    /* NULL algo16 is disallowed */
    CHECK(nonce_function_bip340(nonce, msg, key, pk, NULL, NULL) == 0);
    /* Empty algo16 is fine */
    memset(algo16, 0x00, 16);
    CHECK(nonce_function_bip340(nonce, msg, key, pk, algo16, NULL) == 1);
    /* algo16 with terminating null bytes is fine */
    algo16[1] = 65;
    CHECK(nonce_function_bip340(nonce, msg, key, pk, algo16, NULL) == 1);
    /* Other algo16 is fine */
    memset(algo16, 0xFF, 16);
    CHECK(nonce_function_bip340(nonce, msg, key, pk, algo16, NULL) == 1);

    /* NULL aux_rand argument is allowed. */
    CHECK(nonce_function_bip340(nonce, msg, key, pk, algo16, NULL) == 1);
}

void run_schnorrsig_tests(void) {
    run_nonce_function_bip340_tests();
}

#endif
