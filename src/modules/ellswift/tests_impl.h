/***********************************************************************
 * Copyright (c) 2022 Pieter Wuile                                     *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 ***********************************************************************/

#ifndef SECP256K1_MODULE_ELLSWIFT_TESTS_H
#define SECP256K1_MODULE_ELLSWIFT_TESTS_H

#include "../../../include/secp256k1_ellswift.h"

void run_ellswift_tests(void) {
    int i = 0;
    /* Verify that secp256k1_ellswift_encode + decode roundtrips. */
    for (i = 0; i < 1000 * count; i++) {
        unsigned char rnd32[32];
        unsigned char ell64[64];
        secp256k1_ge g, g2;
        secp256k1_xonly_pubkey pubkey, pubkey2;
        random_group_element_test(&g);
        secp256k1_xonly_pubkey_save(&pubkey, &g);
        secp256k1_testrand256(rnd32);
        secp256k1_ellswift_encode(ctx, ell64, rnd32, &pubkey);
        secp256k1_ellswift_decode(ctx, &pubkey2, ell64);
        secp256k1_xonly_pubkey_load(ctx, &g2, &pubkey2);
        CHECK(check_fe_equal(&g.x, &g2.x));
    }
}

#endif
