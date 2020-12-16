/***********************************************************************
 * Copyright (c) 2020 Peter Dettman, Pieter Wuille                     *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K1_MODINV64_IMPL_H
#define SECP256K1_MODINV64_IMPL_H

#include "modinv64.h"

#include "util.h"

/* This file implements modular inversion based on the paper "Fast constant-time gcd computation and
 * modular inversion" by Daniel J. Bernstein and Bo-Yin Yang.
 *
 * See https://gcd.cr.yp.to/papers.html#safegcd for the paper. The references below are for the Date:
 * 2019.04.13 version.
 *
 * Below is an explanation of the implementation, building up the algorithm in Python3 step by step.
 *
 *
 * 1. Computing the Greatest Common Divisor (GCD) using divsteps
 * -------------------------------------------------------------
 *
 * The algorithm from the paper, at a very high level, is this:
 *
 * def gcd(f, g):
 *     """Compute the GCD of an odd integer f and another integer g."""
 *     assert f & 1  # require f to be odd
 *     delta = 1     # additional state variable
 *     while g != 0:
 *         assert f & 1  # f will be odd in every iteration
 *         if delta > 0 and g & 1:
 *             delta, f, g = 1 - delta, g, (g - f) // 2
 *         elif g & 1:
 *             delta, f, g = 1 + delta, f, (g + f) // 2
 *         else:
 *             delta, f, g = 1 + delta, f, (g    ) // 2
 *     return abs(f)
 *
 * It computes the greatest common divisor of an odd integer f and any integer g. Its inner loop
 * keeps rewriting the variables f and g alongside a state variable delta that starts at 1, until
 * g=0 is reached. At that point, |f| gives the GCD. Each of the transitions in the loop is called a
 * "division step" (referred to as divstep in what follows).
 *
 * For example, gcd(21, 14) would be computed as:
 * - Start with delta=1 f=21 g=14
 * - Take the third branch: delta=2 f=21 g=7
 * - Take the first branch: delta=-1 f=7 g=-7
 * - Take the second branch: delta=0 f=7 g=0
 * - The answer |f| = 7.
 *
 * Why it works:
 * - Divsteps can be decomposed into two steps (see paragraph 8.2 in the paper):
 *   - (a) If g is odd, replace (f,g) with (g,g-f) or (f,g+f), resulting in an even g.
 *   - (b) Replace (f,g) with (f,g//2) (where g is guaranteed to be even).
 * - Neither of those two operations change the GCD:
 *   - For (a), assume gcd(f,g)=c, then it must be the case that f=a*c and g=b*c for some integers a
 *     and b. As (g,g-f)=(b*c,(b-a)*c) and (f,f+g)=(a*c,(a+b)*c), the result clearly still has
 *     common factor c. Reasoning in the other direction shows that no common factor can be added by
 *     doing so either.
 *   - For (b), we know that f is odd, so gcd(f,g) clearly has no factor 2, and we can remove
 *     it from g.
 * - The algorithm will eventually converge to g=0. This is proven in the paper (see theorem G.3).
 * - It follows that eventually we find a final value f' for which gcd(f,g) = gcd(f',0). As the
 *   gcd of f' and 0 is |f'| by definition, that is our answer.
 *
 * Compared to more traditional GCD algorithms, this one has the property of only ever looking at
 * the low-order bits of the variables to decide the next steps, and being easy to make
 * constant-time (in more low-level languages than Python). The delta parameter is necessary to
 * guide the algorithm towards shrinking the numbers' magnitudes without explicitly needing to look
 * at high order bits.
 *
 * Properties that will become important later:
 * - Performing more divsteps than needed is not a problem, as f does not change anymore after g=0.
 * - Only even numbers are divided by 2. This means that when reasoning about it algebraically we
 *   do not need to worry about rounding.
 * - At every point during the algorithm's execution the next N steps only depend on the bottom N
 *   bits of f and g, and on delta.
 *
 *
 * 2. From GCDs to modular inverses
 * --------------------------------
 *
 * We want an algorithm to compute the inverse a of x modulo M, i.e. the number a such that a*x=1
 * mod M. This inverse only exists if the GCD of x and M is 1, but that is always the case if M is
 * prime and 0 < x < M. In what follows, assume that the modular inverse exists. To find that
 * inverse, it turns out this can be computed as a side effect of computing the GCD by keeping track
 * of how the internal variables can be written as linear combinations of the inputs at every step.
 * Since the GCD is 1, such an algorithm will compute numbers a and b such that a*x + b*M = 1.
 * Taking that expression mod M gives a*x mod M = 1, and we see that a is the modular inverse of x
 * mod M.
 *
 * A similar approach can be used to calculate modular inverses using the divsteps-based GCD
 * algorithm shown above, if the modulus M is odd. To do so, compute gcd(f=M,g=x), while keeping
 * track of extra variables d and e, for which at every step d = f/x (mod M) and e = g/x (mod M).
 * f/x here means the number which multiplied with x gives f mod M. As f and g are initialized to M
 * and x respectively, d and e just start off being 0 (M/x mod M = 0/x mod M = 0) and 1 (x/x mod M
 * = 1).
 *
 * def div2(M, x):
 *     """Helper routine to compute x/2 mod M (where M is odd)."""
 *     assert M & 1
 *     if x & 1: # If x is odd, make it even by adding M.
 *         x += M
 *     # x must be even now, so a clean division by 2 is possible.
 *     return x // 2
 *
 * def modinv(M, x):
 *     """Compute the inverse of x mod M (given that it exists, and M is odd)."""
 *     assert M & 1
 *     delta, f, g, d, e = 1, M, x, 0, 1
 *     while g != 0:
 *         # Note that while division by two for f and g is only ever done on even inputs, this is
 *         # not true for d and e, so we need the div2 helper function.
 *         if delta > 0 and g & 1:
 *             delta, f, g, d, e = 1 - delta, g, (g - f) // 2, e, div2(M, e - d)
 *         elif g & 1:
 *             delta, f, g, d, e = 1 + delta, f, (g + f) // 2, d, div2(M, e + d)
 *         else:
 *             delta, f, g, d, e = 1 + delta, f, (g    ) // 2, d, div2(M, e    )
 *         # Verify that the invariants d=f/x mod M, e=g/x mod M are maintained.
 *         assert f % M == (d * x) % M
 *         assert g % M == (e * x) % M
 *     assert f == 1 or f == -1  # |f| is the GCD, it must be 1
 *     # Because of invariant d = f/x (mod M), 1/x = d/f (mod M). As |f|=1, d/f = d*f.
 *     return (d * f) % M
 *
 * Also note that this approach to track d and e throughout the computation to determine the inverse
 * is different from the paper. There (see paragraph 12.1 in the paper) a transition matrix for the
 * entire computation is determined (see section 3 below) and the inverse is computed from that.
 * The approach here avoids the need for 2x2 matrix multiplications of various sizes, and appears to
 * be faster at the level of optimization we're able to do in C.
 *
 *
 * 3. Batching multiple divsteps
 * -----------------------------
 *
 * Every divstep can be expressed as a matrix multiplication, applying a transition matrix (1/2*t)
 * to both vectors [f, g] and [d, e] (see paragraph 8.1 in the paper):
 *
 *   t = [ u,  v ]
 *       [ q,  r ]
 *
 *   [ out_f ] = (1/2 * t) * [ in_f ]
 *   [ out_g ] =             [ in_g ]
 *
 *   [ out_d ] = (1/2 * t) * [ in_d ]  (mod M)
 *   [ out_e ]               [ in_e ]
 *
 * where (u, v, q, r) is (0, 2, -1, 1), (2, 0, 1, 1), or (2, 0, 0, 1), depending on which branch is
 * taken. As above, the resulting f and g are always integers.
 *
 * Performing multiple divsteps corresponds to a multiplication with the product of all the
 * individual divsteps' transition matrices. As each transition matrix consists of integers
 * divided by 2, the product of these matrices will consist of integers divided by 2^N (see also
 * theorem 9.2 in the paper). These divisions are expensive when updating d and e, so we delay
 * them: we compute the integer coefficients of the combined transition matrix scaled by 2^N, and
 * do one division by 2^N as a final step:
 *
 * def divsteps_n_matrix(delta, f, g):
 *     """Compute delta and transition matrix t after N divsteps (multiplied by 2^N)."""
 *     u, v, q, r = 1, 0, 0, 1 # start with identity matrix
 *     for _ in range(N):
 *         if delta > 0 and g & 1:
 *             delta, f, g, u, v, q, r = 1 - delta, g, (g - f) // 2, 2*q, 2*r, q-u, r-v
 *         elif g & 1:
 *             delta, f, g, u, v, q, r = 1 + delta, f, (g + f) // 2, 2*u, 2*v, q+u, r+v
 *         else:
 *             delta, f, g, u, v, q, r = 1 + delta, f, (g    ) // 2, 2*u, 2*v, q  , r
 *     return delta, (u, v, q, r)
 *
 * As the branches in the divsteps are completely determined by the bottom N bits of f and g, this
 * function to compute the transition matrix only needs to see those bottom bits. Furthermore all
 * intermediate results and outputs fit in (N+1)-bit numbers (unsigned for f and g; signed for u, v,
 * q, and r) (see also paragraph 8.3 in the paper). This means that an implementation using 64-bit
 * integers could set N=62 and compute the full transition matrix for 62 steps at once without any
 * big integer arithmetic at all. This is the reason why this algorithm is efficient: it only needs
 * to update the full-size f, g, d, and e numbers once every N steps.
 *
 * We still need functions to compute:
 *
 *   [ out_f ] = (1/2^N * [ u,  v ]) * [ in_f ]
 *   [ out_g ]   (        [ q,  r ])   [ in_g ]
 *
 *   [ out_d ] = (1/2^N * [ u,  v ]) * [ in_d ]  (mod M)
 *   [ out_e ]   (        [ q,  r ])   [ in_e ]
 *
 * For f and g that's easy:
 *
 * def update_fg(f, g, t):
 *     """Multiply matrix t/2^N with [f, g]."""
 *     u, v, q, r = t
 *     cf, cg = u*f + v*g, q*f + r*g
 *     # (t / 2^N) should cleanly apply to [f,g] so the result of t*[f,g] should have N zero
 *     # bottom bits.
 *     assert cf % 2**N == 0
 *     assert cg % 2**N == 0
 *     return cf >> N, cg >> N
 *
 * To do the same for d and e, we need an equivalent of the div2 function for division by 2^N mod M.
 * This is easy if we have precomputed M^-1 mod 2^N (which always exists for odd M):
 *
 * def div2n(M, Mi, x):
 *     """Compute x/2^N mod M, given Mi = 1/M mod 2^N."""
 *     assert (M * Mi) % 2**N == 1
 *     # Find a factor m such that m*M has the same bottom N bits as x. We want:
 *     #     (m * M) mod 2^N = x mod 2^N
 *     # <=> m mod 2^N = (x / M) mod 2^N
 *     # <=> m mod 2^N = (x * Mi) mod 2^N
 *     m = (Mi * x) % 2**N
 *     # Subtract that multiple from x, cancelling its bottom N bits.
 *     x -= m * M
 *     # Now a clean division by 2^N is possible.
 *     assert x % 2**N == 0
 *     return (x >> N) % M
 *
 * def update_de(d, e, t, M, Mi):
 *     """Multiply matrix t/2^N with [d, e], modulo M."""
 *     u, v, q, r = t
 *     cd, ce = u*d + v*e, q*d + r*e
 *     return div2n(M, Mi, cd), div2n(M, Mi, ce)
 *
 * With all of those, we can write a version of modinv that performs N divsteps at once:
 *
 * def modinv(M, Mi, x):
 *     """Compute the modular inverse of x mod M, given Mi=1/M mod 2^N."""
 *     assert M & 1
 *     delta, f, g, d, e = 1, M, x, 0, 1
 *     while g != 0:
 *         # Compute the delta and transition matrix t for the next N divsteps (this only needs
 *         # (N+1)-bit signed integer arithmetic).
 *         delta, t = divsteps_n_matrix(delta, f % 2**N, g % 2**N)
 *         # Apply the transition matrix t to [f, g]:
 *         f, g = update_fg(f, g, t)
 *         # Apply the transition matrix t to [d, e]:
 *         d, e = update_de(d, e, t, M, Mi)
 *     return (d * f) % M
 *
 * This means that in practice we'll always perform a multiple of N divsteps. This is not a problem
 * because once g=0, further divsteps do not affect f, g, d, or e anymore (only delta keeps
 * increasing). For variable time code such excess iterations will be mostly optimized away in
 * section 6.
 *
 *
 * 4. Avoiding modulus operations
 * ------------------------------
 *
 * So far, there are two places where we compute a remainder of big numbers modulo M: at the end of
 * div2n in every update_de, and at the very end of modinv after potentially negating d due to the
 * sign of f. These are relatively expensive operations when done generically.
 *
 * To deal with the modulus operation in div2n, we simply stop requiring d and e to be in range
 * [0,M) all the time. Let's start by inlining div2n into update_de, and dropping the modulus
 * operation at the end:
 *
 * def update_de(d, e, t, M, Mi):
 *     """Multiply matrix t/2^N with [d, e] mod M, given Mi=1/M mod 2^N."""
 *     u, v, q, r = t
 *     cd, ce = u*d + v*e, q*d + r*e
 *     # Cancel out bottom N bits of cd and ce.
 *     md = -((Mi * cd) % 2**N)
 *     me = -((Mi * ce) % 2**N)
 *     cd += md * M
 *     ce += me * M
 *     # And cleanly divide by 2**N.
 *     return cd >> N, ce >> N
 *
 * Let's look at bounds on the ranges of these numbers. It can be shown that |u|+|v| and |q|+|r|
 * never exceed 2^N (see paragraph 8.3 in the paper), and thus a multiplication with t will have
 * outputs whose absolute values are at most 2^N times the maximum absolute input value. In case the
 * inputs d and e are in (-M,M), which is certainly true for the initial values d=0 and e=1 assuming
 * M > 1, the multiplication results in numbers in range (-2^N*M,2^N*M). Subtracting up to 2^N-1
 * times M to cancel out N bits brings that up to slightly less than (-2^(N+1)*M,2^N*M), and
 * dividing by 2^N at the end takes it to (-2*M,M). Another application of update_de would take that
 * to (-3*M,2*M), and so forth. This progressive expansion of the variables' ranges can be
 * counteracted by incrementing d and e by M whenever they're negative:
 *
 *     ...
 *     if d < 0:
 *         d += M
 *     if e < 0:
 *         e += M
 *     cd, ce = u*d + v*e, q*d + r*e
 *     # Cancel out bottom N bits of cd and ce.
 *     ...
 *
 * With inputs in (-2*M,M), they will first be shifted into range (-M,M), which means that the
 * output will again be in (-2*M,M), and this remains the case regardless of how many update_de
 * invocations there are. In what follows, we will try to make this more efficient.
 *
 * Note that increasing d by M is equal to incrementing cd by u*M and ce by q*M. Similarly,
 * increasing e by M is equal to incrementing cd by v*M and ce by r*M. So we could instead write:
 *
 *     ...
 *     cd, ce = u*d + v*e, q*d + r*e
 *     # Perform the equivalent of incrementing d, e by M when they're negative.
 *     if d < 0:
 *         cd += u*M
 *         ce += q*M
 *     if e < 0:
 *         cd += v*M
 *         ce += r*M
 *     # Cancel out bottom N bits of cd and ce.
 *     md = -((Mi * cd) % 2**N)
 *     me = -((Mi * ce) % 2**N)
 *     cd += md * M
 *     ce += me * M
 *     ...
 *
 * Now note that we have two steps of corrections to cd and ce that add multiples of M: this
 * increment, and the decrement that cancels out bottom bits. The second one depends on the first
 * one, but they can still be efficiently combined by only computing the bottom bits of cd and ce
 * at first, and using that to compute the final md, me values:
 *
 * def update_de(d, e, t, M, Mi):
 *     """Multiply matrix t/2^N with [d, e], modulo M."""
 *     u, v, q, r = t
 *     md, me = 0, 0
 *     # Compute what multiples of M to add to cd and ce.
 *     if d < 0:
 *         md += u
 *         me += q
 *     if e < 0:
 *         md += v
 *         me += r
 *     # Compute bottom N bits of t*[d,e] + M*[md,me].
 *     cd, ce = (u*d + v*e + md*M) % 2**N, (q*d + r*e + me*M) % 2**N
 *     # Correct md and me such that the bottom N bits of t*[d,e] + M*[md,me] are zero.
 *     md -= (Mi * cd) % 2**N
 *     me -= (Mi * ce) % 2**N
 *     # Do the full computation.
 *     cd, ce = u*d + v*e + md*M, q*d + r*e + me*M
 *     # And cleanly divide by 2**N.
 *     return cd >> N, ce >> N
 *
 * One last optimization: we can avoid the md*M and me*M multiplications in the bottom bits of cd
 * and ce by moving them to the md and me correction:
 *
 *     ...
 *     # Compute bottom N bits of t*[d,e].
 *     cd, ce = (u*d + v*e) % 2**N, (q*d + r*e) % 2**N
 *     # Correct md and me such that the bottom N bits of t*[d,e]+M*[md,me] are zero.
 *     # Note that this is not the same as {md = (Mi * cd) % 2**N} etc. That would also result in N
 *     # zero bottom bits, but isn't guaranteed to be a reduction of [0,2^N) compared to the
 *     # previous md and me values, and thus would violate our bounds analysis.
 *     md -= (Mi*cd + md) % 2**N
 *     me -= (Mi*ce + me) % 2**N
 *     ...
 *
 * The resulting function takes d and e in range (-2*M,M) as inputs, and outputs values in the same
 * range. That also means that the d value at the end of modinv will be in that range, while we want
 * a result in [0,M). To do that, we need a normalization function. It's easy to integrate the
 * conditional negation of d (based on the sign of f) into it as well:
 *
 * def normalize(sign, v, M):
 *     """Compute sign*v mod M, where v is in range (-2*M,M); output in [0,M)."""
 *     assert sign == 1 or sign == -1
 *     # v in (-2*M,M)
 *     if v < 0:
 *         v += M
 *     # v in (-M,M). Now multiply v with sign (which can only be 1 or -1).
 *     if sign == -1:
 *         v = -v
 *     # v in (-M,M)
 *     if v < 0:
 *         v += M
 *     # v in [0,M)
 *     return v
 *
 * And calling it in modinv is simply:
 *
 *    ...
 *    return normalize(f, d, M)
 *
 *
 * 5. Constant-time operation
 * --------------------------
 *
 * The primary selling point of the algorithm is fast constant-time operation. What code flow still
 * depends on the input data so far?
 *
 * - The number of iterations of the while g != 0 loop in modinv.
 * - The branches inside divsteps_n_matrix.
 * - The sign checks in update_de
 * - The sign checks in normalize
 *
 * To make the while loop in modinv constant time it can be replaced with a constant number of
 * iterations. The paper proves (Theorem 11.2) that 741 divsteps are sufficient for any 256-bit
 * inputs, and https://github.com/sipa/safegcd-bounds shows that the slightly better bound 724 is
 * sufficient even. Given that every loop iteration performs N divsteps, it will run a total of
 * ceil(724/N) times.
 *
 * To deal with the branches in divsteps_n_matrix we will replace them with constant-time bitwise
 * operations (and hope the C compiler isn't smart enough to turn them back into branches; see
 * valgrind_ctime_test.c for automated tests that this isn't the case). To do so, observe that a
 * divstep can be written instead as (compare to the inner loop of gcd() in section 1).
 *
 *     x = -f if delta > 0 else f         # set x equal to (input) -f or f
 *     if g & 1:
 *         g += x                         # set g to (input) g-f or g+f
 *         if delta > 0:
 *             delta = -delta
 *             f += g                     # set f to (input) g (note that g was set to g-f before)
 *     delta += 1
 *     g >>= 1
 *
 * To convert the above to bitwise operations, we rely on a trick to negate conditionally: per the
 * definition of negative numbers in two's complement, (-v == ~v + 1) holds for every number v. As
 * -1 in two's complement is all 1 bits, bitflipping can be expressed as xor with -1. It follows
 * that (-v == (v ^ -1) - (-1)). Thus, if we have a variable c that takes on values 0 or -1, then
 * ((v ^ c) - c) is v if c=0 and -v if c=-1.
 *
 * Using this we can write:
 *     x = -f if delta > 0 else f
 * in constant-time form as:
 *     c1 = (-delta) >> 63
 *     # Conditionally negate f based on c1:
 *     x = (f ^ c1) - c1
 *
 * To use that trick, we need a helper mask variable c1 that resolves the condition delta>0 to -1
 * (if true) or 0 (if false). We compute c1 using right shifting, which always rounds down (in
 * Python, and also in C under the assumption of a typical two's complement system; see
 * assumptions.h for tests that this is the case).
 *
 * Using the facts that x&0=0 and x&(-1)=x (on two's complement systems again), we can write:
 *     if g & 1:
 *         g += x
 * as:
 *     # Compute c2=0 if g is even and c2=-1 if g is odd.
 *     c2 = -(g & 1)
 *     # This masks out x if g is even, and leaves x be if g is odd.
 *     g += x & c2
 *
 * Using the conditional negation trick again we can write:
 *     if g & 1:
 *         if delta > 0:
 *             delta = -delta
 * as:
 *     # Compute c3=-1 if g is odd and delta>0, and 0 otherwise.
 *     c3 = c1 & c2
 *     # Conditionally negate delta based on c3:
 *     delta = (delta ^ c3) - c3
 *
 * Finally:
 *     if g & 1:
 *         if delta > 0:
 *             f += g
 * becomes:
 *     f += g & c3
 *
 * It turns out that this can be implemented more efficiently by applying the substitution
 * eta=-delta. In this representation, negating delta corresponds to negating eta, and incrementing
 * delta corresponds to decrementing eta. This allows us to remove the negation in the c1
 * computation:
 *
 *     # Compute a mask c1 for eta < 0, and compute the conditional negation x of f:
 *     c1 = eta >> 63
 *     x = (f ^ c1) - c1
 *     # Compute a mask c2 for odd g, and conditionally add x to g:
 *     c2 = -(g & 1)
 *     g += x & c2
 *     # Compute a mask c for (eta < 0) and odd (input) g, and use it to conditionally negate eta,
 *     # and add g to f:
 *     c3 = c1 & c2
 *     eta = (eta ^ c3) - c3
 *     f += g & c3
 *     # Incrementing delta corresponds to decrementing eta.
 *     eta -= 1
 *     g >>= 1
 *
 * By replacing the loop in divsteps_n_matrix with a variant of the divstep code above (extended to
 * also apply all f operations to u, v and all g operations to q, r), a constant-time version of
 * divsteps_n_matrix is obtained. The full code will be in section 7.
 *
 * These bit fiddling tricks can also be used to make the conditional negations and additions in
 * update_de and normalize constant-time.
 *
 *
 * 6. Variable-time optimizations
 * ------------------------------
 *
 * In section 5, we modified the divsteps_n_matrix function (and a few others) to be constant time.
 * Constant time operations are only necessary when computing modular inverses of secret data. In
 * other cases, it slows down calculations unnecessarily. In this section, we will construct a
 * faster non-constant time divsteps_n_matrix function.
 *
 * To do so, first consider yet another way of writing the inner loop of divstep operations in
 * gcd() from section 1. This decomposition is also explained in the paper in section 8.2.
 *
 * for _ in range(N):
 *     if g & 1 and eta < 0:
 *         eta, f, g = -eta, g, -f
 *     if g & 1:
 *         g += f
 *     eta -= 1
 *     g >>= 1
 *
 * Whenever g is even, the loop only shifts g down and decreases eta. When g ends in multiple zero
 * bits, these iterations can be consolidated into one step. This requires counting the bottom zero
 * bits efficiently, which is possible on most platforms; it is abstracted here as the function
 * count_trailing_zeros.
 *
 * i = N # divsteps left to do
 * while True:
 *     # Get rid of all bottom zeros at once. In the first iteration, g may be odd and the following
 *     # lines have no effect (until "if eta < 0").
 *     zeros = min(i, count_trailing_zeros(g))
 *     eta -= zeros
 *     g >>= zeros
 *     i -= zeros
 *     if i == 0:
 *         break
 *     # We know g is odd now
 *     if eta < 0:
 *         eta, f, g = -eta, g, -f
 *     g += f
 *     # g is even now, and the eta decrement and g shift will happen in the next loop.
 *
 * We can now remove multiple bottom 0 bits from g at once, but still need a full iteration whenever
 * there is a bottom 1 bit. In what follows, we will get rid of multiple 1 bits simultaneously as
 * well.
 *
 * Observe that as long as eta is >= 0, the loop does not modify f. Instead, it cancels out bottom
 * bits of g and shifts them out, and decreases eta and i accordingly - interrupting only when eta
 * becomes negative, or when i reaches 0. Combined, this is equivalent to adding a multiple of f to
 * g to cancel out multiple bottom bits, and then shifting them out.
 *
 * It is easy to find what that multiple is: we want a number w such that g+w*f has a few bottom
 * zero bits. If that number of bits is L, we want g+w*f mod 2^L = 0, or w = -g/f mod 2^L. Since f
 * is odd, such a w exists for any L. L cannot be more than i steps (as we'd finish the loop before
 * doing more) or more than eta+1 steps (as we'd run {eta, f, g = -eta, g, f} at that point), but
 * apart from that, we're only limited by the complexity of computing w.
 *
 * This code demonstrates how to cancel up to 4 bits per step:
 *
 * NEGINV16 = [15, 5, 3, 9, 7, 13, 11, 1] # NEGINV16[n//2] = (-n)^-1 mod 16, for odd n
 * i = N
 * while True:
 *     zeros = min(i, count_trailing_zeros(g))
 *     eta -= zeros
 *     g >>= zeros
 *     i -= zeros
 *     if i == 0:
 *         break
 *     # We know g is odd now
 *     if eta < 0:
 *         eta, f, g = -eta, g, f
 *     # Compute limit on number of bits to cancel
 *     limit = min(min(eta + 1, i), 4)
 *     # Compute w = -g/f mod 2**limit, using the table value for -1/f mod 2**4. Note that f is
 *     # always odd, so its inverse modulo a power of two always exists.
 *     w = (g * NEGINV16[(f & 15) // 2]) % (2**limit)
 *     # As w = -g/f mod (2**limit), g+w*f mod 2**limit = 0 mod 2**limit.
 *     g += w * f
 *     assert g % (2**limit) == 0
 *     # The next iteration will now shift out at least limit bottom zero bits from g.
 *
 * By using a bigger table more bits can be cancelled at once. The table can also be implemented
 * as a formula:
 *  - Instead of a 3-bit table:
 *    * (-f) or (f ^ 6)
 *  - Instead of a 4-bit table:
 *    * (1 - f * (f + 1))
 *    * (-(f + (((f + 1) & 4) << 1)))
 *  - For larger tables the following technique can be used: if w=-f^-1 mod 2^l, then w*(w*f+2) is
 *    -f^-1 mod 2^(2*l). This allows extending the previous formulas (or tables). In particular we
 *    have this 6-bit function (based on the 3-bit function above):
 *    * (f * (f * f - 2))
 *
 * This loop, again extended to also handle u, v, q, and r alongside f and g, placed in
 * divsteps_n_matrix, gives a significantly faster, but non-constant time version.
 *
 *
 * 7. Final Python version
 * -----------------------
 *
 * All together we need the following functions:
 *
 * - A way to compute the transition matrix in constant time, using the divsteps_n_matrix function
 *   from section 2, but with its loop replaced by a variant of the constant-time divstep from
 *   section 5, extended to handle u, v, q, r:
 *
 * def divsteps_n_matrix(eta, f, g):
 *     """Compute eta and transition matrix t after N divsteps (multiplied by 2^N)."""
 *     u, v, q, r = 1, 0, 0, 1 # start with identity matrix
 *     for _ in range(N):
 *         c1 = eta >> 63
 *         # Compute x, y, z as conditionally-negated versions of f, u, v.
 *         x, y, z = (f ^ c1) - c1, (u ^ c1) - c1, (v ^ c1) - c1
 *         c2 = -(g & 1)
 *         # Conditionally add x, y, z to g, q, r.
 *         g, q, r = g + (x & c2), q + (y & c2), r + (z & c2)
 *         c1 &= c2                     # reusing c1 here for the earlier c3 variable
 *         eta = (eta ^ c1) - (c1 + 1)  # inlining the unconditional eta decrement here
 *         # Conditionally add g, q, r to f, u, v.
 *         f, u, v = f + (g & c1), u + (q & c1), v + (r & c1)
 *         # When shifting g down, don't shift q, r, as we construct a transition matrix multiplied
 *         # by 2^N. Instead, shift f's coefficients u and v up.
 *         g, u, v = g >> 1, u << 1, v << 1
 *     return eta, (u, v, q, r)
 *
 * - The functions to update f and g, and d and e, from section 2 and section 4, with the constant-time
 *   changes to update_de from section 5:
 *
 * def update_fg(f, g, t):
 *     """Multiply matrix t/2^N with [f, g]."""
 *     u, v, q, r = t
 *     cf, cg = u*f + v*g, q*f + r*g
 *     return cf >> N, cg >> N
 *
 * def update_de(d, e, t, M, Mi):
 *     """Multiply matrix t/2^N with [d, e], modulo M."""
 *     u, v, q, r = t
 *     d_sign, e_sign = d >> 257, e >> 257
 *     md, me = (u & d_sign) + (v & e_sign), (q & d_sign) + (r & e_sign)
 *     cd, ce = (u*d + v*e) % 2**N, (q*d + r*e) % 2**N
 *     md -= (Mi*cd + md) % 2**N
 *     me -= (Mi*ce + me) % 2**N
 *     cd, ce = u*d + v*e + Mi*md, q*d + r*e + Mi*me
 *     return cd >> N, ce >> N
 *
 * - The normalize function from section 4, made constant time as well:
 *
 * def normalize(sign, v, M):
 *     """Compute sign*v mod M, where v in (-2*M,M); output in [0,M)."""
 *     v_sign = v >> 257
 *     # Conditionally add M to v.
 *     v += M & v_sign
 *     c = (sign - 1) >> 1
 *     # Conditionally negate v.
 *     v = (v ^ c) - c
 *     v_sign = v >> 257
 *     # Conditionally add M to v again.
 *     v += M & v_sign
 *     return v
 *
 * - And finally the modinv function too, adapted to use eta instead of delta, and using the fixed
 *   iteration count from section 5:
 *
 * def modinv(M, Mi, x):
 *     """Compute the modular inverse of x mod M, given Mi=1/M mod 2^N."""
 *     eta, f, g, d, e = -1, M, x, 0, 1
 *     for _ in range((724 + N - 1) // N):
 *         eta, t = divsteps_n_matrix(-eta, f % 2**N, g % 2**N)
 *         f, g = update_fg(f, g, t)
 *         d, e = update_de(d, e, t, M, Mi)
 *     return normalize(f, d, M)
 *
 * - To get a variable time version, replace the divsteps_n_matrix function with one that uses the
 *   divsteps loop from section 5, and a modinv version that calls it without the fixed iteration
 *   count:
 *
 * NEGINV16 = [15, 5, 3, 9, 7, 13, 11, 1] # NEGINV16[n//2] = (-n)^-1 mod 16, for odd n
 * def divsteps_n_matrix_var(eta, f, g):
 *     """Compute eta and transition matrix t after N divsteps (multiplied by 2^N)."""
 *     u, v, q, r = 1, 0, 0, 1
 *     i = N
 *     while True:
 *         zeros = min(i, count_trailing_zeros(g))
 *         eta, i = eta - zeros, i - zeros
 *         g, u, v = g >> zeros, u << zeros, v << zeros
 *         if i == 0:
 *             break
 *         if eta < 0:
 *             eta, f, u, v, g, q, r = -eta, g, q, r, -f, -u, -v
 *         limit = min(min(eta + 1, i), 4)
 *         w = (g * NEGINV16[(f & 15) // 2]) % (2**limit)
 *         g, q, r = g + w*f, q + w*u, r + w*v
 *     return eta, (u, v, q, r)
 *
 * def modinv_var(M, Mi, x):
 *     """Compute the modular inverse of x mod M, given Mi = 1/M mod 2^N."""
 *     eta, f, g, d, e = -1, M, x, 0, 1
 *     while g != 0:
 *         eta, t = divsteps_n_matrix_var(eta, f % 2**N, g % 2**N)
 *         f, g = update_fg(f, g, t)
 *         d, e = update_de(d, e, t, M, Mi)
 *     return normalize(f, d, Mi)
 *
 *
 * 8. C implementation
 * -------------------
 *
 * What follows is a C implementation of effectively the Python code from section 7, with N=62, and
 * the following changes:
 *
 * - Representing large integers using 5 62-bit signed limbs that callers need to convert their
 *   value from/to. Using 62-bit limbs means shifts by 62 bits are very efficient, and the extra
 *   space allows faster operations by delaying carries/borrows in some cases.
 *
 * - Several modulo operations in the Python code are modulo a power of two. These can be replaced
 *   with a bitwise AND with ((1 << bits) - 1).
 *
 * - Similarly, if an entire expression involving multiplications and additions is computed modulo
 *   a power of two, that means only the bottom bits of the inputs and intermediary results is
 *   needed.
 */

#ifdef VERIFY
/* Helper function to compute the absolute value of an int64_t.
 * (we don't use abs/labs/llabs as it depends on the int sizes). */
static int64_t secp256k1_modinv64_abs(int64_t v) {
    VERIFY_CHECK(v > INT64_MIN);
    if (v < 0) return -v;
    return v;
}

static const secp256k1_modinv64_signed62 SECP256K1_SIGNED62_ONE = {{1}};

/* Compute a*factor and put it in r. All but the top limb in r will be in range [0,2^62). */
static void secp256k1_modinv64_mul_62(secp256k1_modinv64_signed62 *r, const secp256k1_modinv64_signed62 *a, int64_t factor) {
    const int64_t M62 = (int64_t)(UINT64_MAX >> 2);
    int128_t c = 0;
    int i;
    for (i = 0; i < 4; ++i) {
        c += (factor == 1) ? (int128_t)a->v[i] : (int128_t)a->v[i] * factor;
        r->v[i] = (int64_t)c & M62; c >>= 62;
    }
    c += (factor == 1) ? (int128_t)a->v[4] : (int128_t)a->v[4] * factor;
    VERIFY_CHECK(c == (int64_t)c);
    r->v[4] = (int64_t)c;
}

/* Compare af with b*bf. */
static int secp256k1_modinv64_mul_cmp_62(const secp256k1_modinv64_signed62 *a, const secp256k1_modinv64_signed62 *b, int64_t factor) {
    int i;
    secp256k1_modinv64_signed62 am, bm;
    secp256k1_modinv64_mul_62(&am, a, 1);
    secp256k1_modinv64_mul_62(&bm, b, factor);
    for (i = 4; i >= 0; --i) {
        if (i != 4) {
            VERIFY_CHECK(am.v[i] >> 62 == 0);
            VERIFY_CHECK(bm.v[i] >> 62 == 0);
        }
        if (am.v[i] < bm.v[i]) return -1;
        if (am.v[i] > bm.v[i]) return 1;
    }
    return 0;
}
#endif

/* Optionally negate a signed62 number in range (-2*modulus,modulus), and add multiple of modulus to
 * bring it to [0,modulus). The input must have limbs in range (-2^62,2^62). The output will have
 * limbs in range [0,2^62). */
static void secp256k1_modinv64_normalize_62(secp256k1_modinv64_signed62 *r, int64_t cond_negate, const secp256k1_modinv64_modinfo *modinfo) {
    const int64_t M62 = (int64_t)(UINT64_MAX >> 2);
    int64_t r0 = r->v[0], r1 = r->v[1], r2 = r->v[2], r3 = r->v[3], r4 = r->v[4];
    int64_t cond_add;

#ifdef VERIFY
    /* Verify that all limbs are in range (-2^62,2^62). */
    int i;
    for (i = 0; i < 5; ++i) {
        VERIFY_CHECK(r->v[i] >= -M62);
        VERIFY_CHECK(r->v[i] <= M62);
    }
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(r, &modinfo->modulus, -2) > 0); /* r > -2*modulus */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(r, &modinfo->modulus, 1) < 0); /* r < modulus */
#endif

    /* In a first step, add the modulus if the input is negative, and then negate if requested.
     * This brings r from range (-2*modulus,modulus) to range (-modulus,modulus). As all input
     * limbs are in range (-2^62,2^62), this cannot overflow an int64_t. */
    cond_add = r4 >> 63;
    r0 += modinfo->modulus.v[0] & cond_add;
    r1 += modinfo->modulus.v[1] & cond_add;
    r2 += modinfo->modulus.v[2] & cond_add;
    r3 += modinfo->modulus.v[3] & cond_add;
    r4 += modinfo->modulus.v[4] & cond_add;
    r0 = (r0 ^ cond_negate) - cond_negate;
    r1 = (r1 ^ cond_negate) - cond_negate;
    r2 = (r2 ^ cond_negate) - cond_negate;
    r3 = (r3 ^ cond_negate) - cond_negate;
    r4 = (r4 ^ cond_negate) - cond_negate;
    /* Propagate the top bits, to bring limbs back to range (-2^62,2^62). */
    r1 += r0 >> 62; r0 &= M62;
    r2 += r1 >> 62; r1 &= M62;
    r3 += r2 >> 62; r2 &= M62;
    r4 += r3 >> 62; r3 &= M62;

    /* In a second step add the modulus again if the result is still negative, bringing
     * r to range [0,modulus). */
    cond_add = r4 >> 63;
    r0 += modinfo->modulus.v[0] & cond_add;
    r1 += modinfo->modulus.v[1] & cond_add;
    r2 += modinfo->modulus.v[2] & cond_add;
    r3 += modinfo->modulus.v[3] & cond_add;
    r4 += modinfo->modulus.v[4] & cond_add;
    /* And propagate again. */
    r1 += r0 >> 62; r0 &= M62;
    r2 += r1 >> 62; r1 &= M62;
    r3 += r2 >> 62; r2 &= M62;
    r4 += r3 >> 62; r3 &= M62;

    r->v[0] = r0;
    r->v[1] = r1;
    r->v[2] = r2;
    r->v[3] = r3;
    r->v[4] = r4;

#ifdef VERIFY
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(r, &modinfo->modulus, 0) >= 0); /* r >= 0 */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(r, &modinfo->modulus, 1) < 0); /* r < P */
    VERIFY_CHECK(r->v[0] >> 62 == 0);
    VERIFY_CHECK(r->v[1] >> 62 == 0);
    VERIFY_CHECK(r->v[2] >> 62 == 0);
    VERIFY_CHECK(r->v[3] >> 62 == 0);
    VERIFY_CHECK(r->v[4] >>  8 == 0);
#endif
}

/* Data type for transition matrices (see section 3 of explanation).
 *
 * t = [ u  v ]
 *     [ q  r ]
 */
typedef struct {
    int64_t u, v, q, r;
} secp256k1_modinv64_trans2x2;

/* Compute the transition matrix and eta for 62 divsteps.
 *
 * Input:  eta: initial eta
 *         f0:  bottom limb of initial f
 *         g0:  bottom limb of initial g
 * Output: t: transition matrix
 * Return: final eta
 *
 * Implements the divsteps_n_matrix function from the explanation.
 */
static uint64_t secp256k1_modinv64_divsteps_62(uint64_t eta, uint64_t f0, uint64_t g0, secp256k1_modinv64_trans2x2 *t) {
    uint64_t u = 1, v = 0, q = 0, r = 1; /* start with identity matrix */
    uint64_t c1, c2, f = f0, g = g0, x, y, z;
    int i;

    for (i = 0; i < 62; ++i) {
        VERIFY_CHECK((f & 1) == 1); /* f must always be odd */
        VERIFY_CHECK((u * f0 + v * g0) == f << i);
        VERIFY_CHECK((q * f0 + r * g0) == g << i);
        /* Compute conditional masks for (eta < 0) and for (g & 1). */
        c1 = (int64_t)eta >> 63;
        c2 = -(g & 1);
        /* Compute x,y,z, conditionally negated versions of f,u,v. */
        x = (f ^ c1) - c1;
        y = (u ^ c1) - c1;
        z = (v ^ c1) - c1;
        /* Conditionally add x,y,z to g,q,r. */
        g += x & c2;
        q += y & c2;
        r += z & c2;
        /* In what follows, c1 is a condition mask for (eta < 0) and (g & 1). */
        c1 &= c2;
        /* Conditionally negate eta, and unconditionally subtract 1. */
        eta = (eta ^ c1) - (c1 + 1);
        /* Conditionally add g,q,r to f,u,v. */
        f += g & c1;
        u += q & c1;
        v += r & c1;
        /* Shifts */
        g >>= 1;
        u <<= 1;
        v <<= 1;
    }
    /* Return data in t and return value. */
    t->u = (int64_t)u;
    t->v = (int64_t)v;
    t->q = (int64_t)q;
    t->r = (int64_t)r;
    return eta;
}

/* Compute the transition matrix and eta for 62 divsteps (variable time).
 *
 * Input:  eta: initial eta
 *         f0:  bottom limb of initial f
 *         g0:  bottom limb of initial g
 * Output: t: transition matrix
 * Return: final eta
 *
 * Implements the divsteps_n_matrix_var function from the explanation.
 */
static uint64_t secp256k1_modinv64_divsteps_62_var(uint64_t eta, uint64_t f0, uint64_t g0, secp256k1_modinv64_trans2x2 *t) {
    uint64_t u = 1, v = 0, q = 0, r = 1; /* Start with identity matrix */
    uint64_t f = f0, g = g0, m, w, x, y, z;
    int i = 62, limit, zeros;

    for (;;) {
        /* Use a sentinel bit to count zeros only up to i. */
        zeros = secp256k1_ctz64_var(g | (UINT64_MAX << i));
        /* Perform zeros divsteps at once; they all just divide g by two. */
        g >>= zeros;
        u <<= zeros;
        v <<= zeros;
        eta -= zeros;
        i -= zeros;
        /* We're done once we've done 62 divsteps. */
        if (i == 0) break;
        VERIFY_CHECK((f & 1) == 1);
        VERIFY_CHECK((g & 1) == 1);
        VERIFY_CHECK((u * f0 + v * g0) == f << (62 - i));
        VERIFY_CHECK((q * f0 + r * g0) == g << (62 - i));
        /* If eta is negative, negate it and replace f,g with g,-f. */
        if ((int64_t)eta < 0) {
            eta = -eta;
            x = f; f = g; g = -x;
            y = u; u = q; q = -y;
            z = v; v = r; r = -z;
            /* Use a formula to cancel out up to 6 bits of g. Also, no more than i can be cancelled
             * out (as we'd be done before that point), and no more than eta+1 can be done as its
             * will flip again once that happens. */
            limit = ((int)eta + 1) > i ? i : ((int)eta + 1);
            /* m is a mask for the bottom min(limit, 6) bits. */
            m = (UINT64_MAX >> (64 - limit)) & 63U;
            /* Find what multiple of f must be added to g to cancel its bottom min(limit, 6)
             * bits. */
            w = (f * g * (f * f - 2)) & m;
        } else {
            /* In this branch, use a simpler formula that only lets us cancel up to 4 bits of g, as
             * eta tends to be smaller here. */
            limit = ((int)eta + 1) > i ? i : ((int)eta + 1);
            /* m is a mask for the bottom min(limit, 4) bits. */
            m = (UINT64_MAX >> (64 - limit)) & 15U;
            /* Find what multiple of f must be added to g to cancel its bottom min(limit, 6)
             * bits. */
            w = f + (((f + 1) & 4) << 1);
            w = (-w * g) & m;
        }
        /* Cancel those bits of g. */
        g += f * w;
        q += u * w;
        r += v * w;
        VERIFY_CHECK((g & m) == 0);
    }
    /* Return data in t and return value. */
    t->u = (int64_t)u;
    t->v = (int64_t)v;
    t->q = (int64_t)q;
    t->r = (int64_t)r;
    return eta;
}

/* Compute (t/2^62) * [d, e] mod modulus, where t is a transition matrix for 62 divsteps.
 *
 * On input and output, d and e are in range (-2*modulus,modulus). All output limbs will be in range
 * (-2^62,2^62).
 *
 * This implements the update_de function from the explanation.
 */
static void secp256k1_modinv64_update_de_62(secp256k1_modinv64_signed62 *d, secp256k1_modinv64_signed62 *e, const secp256k1_modinv64_trans2x2 *t, const secp256k1_modinv64_modinfo* modinfo) {
    const int64_t M62 = (int64_t)(UINT64_MAX >> 2);
    const int64_t d0 = d->v[0], d1 = d->v[1], d2 = d->v[2], d3 = d->v[3], d4 = d->v[4];
    const int64_t e0 = e->v[0], e1 = e->v[1], e2 = e->v[2], e3 = e->v[3], e4 = e->v[4];
    const int64_t u = t->u, v = t->v, q = t->q, r = t->r;
    int64_t md, me, sd, se;
    int128_t cd, ce;
#ifdef VERIFY
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(d, &modinfo->modulus, -2) > 0); /* d > -2*modulus */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(d, &modinfo->modulus, 1) < 0);  /* d <    modulus */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(e, &modinfo->modulus, -2) > 0); /* e > -2*modulus */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(e, &modinfo->modulus, 1) < 0);  /* e <    modulus */
    VERIFY_CHECK((secp256k1_modinv64_abs(u) + secp256k1_modinv64_abs(v)) >= 0); /* |u|+|v| doesn't overflow */
    VERIFY_CHECK((secp256k1_modinv64_abs(q) + secp256k1_modinv64_abs(r)) >= 0); /* |q|+|r| doesn't overflow */
    VERIFY_CHECK((secp256k1_modinv64_abs(u) + secp256k1_modinv64_abs(v)) <= M62 + 1); /* |u|+|v| <= 2^62 */
    VERIFY_CHECK((secp256k1_modinv64_abs(q) + secp256k1_modinv64_abs(r)) <= M62 + 1); /* |q|+|r| <= 2^62 */
#endif
    /* [md,me] start as zero; plus [u,q] if d is negative; plus [v,r] if e is negative. */
    sd = d4 >> 63;
    se = e4 >> 63;
    md = (u & sd) + (v & se);
    me = (q & sd) + (r & se);
    /* Begin computing t*[d,e]. */
    cd = (int128_t)u * d0 + (int128_t)v * e0;
    ce = (int128_t)q * d0 + (int128_t)r * e0;
    /* Correct md,me so that t*[d,e]+modulus*[md,me] has 62 zero bottom bits. */
    md -= (modinfo->modulus_inv62 * (uint64_t)cd + md) & M62;
    me -= (modinfo->modulus_inv62 * (uint64_t)ce + me) & M62;
    /* Update the beginning of computation for t*[d,e]+modulus*[md,me] now md,me are known. */
    cd += (int128_t)modinfo->modulus.v[0] * md;
    ce += (int128_t)modinfo->modulus.v[0] * me;
    /* Verify that the low 62 bits of the computation are indeed zero, and then throw them away. */
    VERIFY_CHECK(((int64_t)cd & M62) == 0); cd >>= 62;
    VERIFY_CHECK(((int64_t)ce & M62) == 0); ce >>= 62;
    /* Compute limb 1 of t*[d,e]+modulus*[md,me], and store it as output limb 0 (= down shift). */
    cd += (int128_t)u * d1 + (int128_t)v * e1;
    ce += (int128_t)q * d1 + (int128_t)r * e1;
    if (modinfo->modulus.v[1]) { /* Optimize for the case where limb of modulus is zero. */
        cd += (int128_t)modinfo->modulus.v[1] * md;
        ce += (int128_t)modinfo->modulus.v[1] * me;
    }
    d->v[0] = (int64_t)cd & M62; cd >>= 62;
    e->v[0] = (int64_t)ce & M62; ce >>= 62;
    /* Compute limb 2 of t*[d,e]+modulus*[md,me], and store it as output limb 1. */
    cd += (int128_t)u * d2 + (int128_t)v * e2;
    ce += (int128_t)q * d2 + (int128_t)r * e2;
    if (modinfo->modulus.v[2]) { /* Optimize for the case where limb of modulus is zero. */
        cd += (int128_t)modinfo->modulus.v[2] * md;
        ce += (int128_t)modinfo->modulus.v[2] * me;
    }
    d->v[1] = (int64_t)cd & M62; cd >>= 62;
    e->v[1] = (int64_t)ce & M62; ce >>= 62;
    /* Compute limb 3 of t*[d,e]+modulus*[md,me], and store it as output limb 2. */
    cd += (int128_t)u * d3 + (int128_t)v * e3;
    ce += (int128_t)q * d3 + (int128_t)r * e3;
    if (modinfo->modulus.v[3]) { /* Optimize for the case where limb of modulus is zero. */
        cd += (int128_t)modinfo->modulus.v[3] * md;
        ce += (int128_t)modinfo->modulus.v[3] * me;
    }
    d->v[2] = (int64_t)cd & M62; cd >>= 62;
    e->v[2] = (int64_t)ce & M62; ce >>= 62;
    /* Compute limb 4 of t*[d,e]+modulus*[md,me], and store it as output limb 3. */
    cd += (int128_t)u * d4 + (int128_t)v * e4;
    ce += (int128_t)q * d4 + (int128_t)r * e4;
    cd += (int128_t)modinfo->modulus.v[4] * md;
    ce += (int128_t)modinfo->modulus.v[4] * me;
    d->v[3] = (int64_t)cd & M62; cd >>= 62;
    e->v[3] = (int64_t)ce & M62; ce >>= 62;
    /* What remains is limb 5 of t*[d,e]+modulus*[md,me]; store it as output limb 4. */
    d->v[4] = (int64_t)cd;
    e->v[4] = (int64_t)ce;
#ifdef VERIFY
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(d, &modinfo->modulus, -2) > 0); /* d > -2*modulus */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(d, &modinfo->modulus, 1) < 0);  /* d <    modulus */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(e, &modinfo->modulus, -2) > 0); /* e > -2*modulus */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(e, &modinfo->modulus, 1) < 0);  /* e <    modulus */
#endif
}

/* Compute (t/2^62) * [f, g], where t is a transition matrix for 62 divsteps.
 *
 * This implements the update_fg function from the explanation.
 */
static void secp256k1_modinv64_update_fg_62(secp256k1_modinv64_signed62 *f, secp256k1_modinv64_signed62 *g, const secp256k1_modinv64_trans2x2 *t) {
    const int64_t M62 = (int64_t)(UINT64_MAX >> 2);
    const int64_t f0 = f->v[0], f1 = f->v[1], f2 = f->v[2], f3 = f->v[3], f4 = f->v[4];
    const int64_t g0 = g->v[0], g1 = g->v[1], g2 = g->v[2], g3 = g->v[3], g4 = g->v[4];
    const int64_t u = t->u, v = t->v, q = t->q, r = t->r;
    int128_t cf, cg;
    /* Start computing t*[f,g]. */
    cf = (int128_t)u * f0 + (int128_t)v * g0;
    cg = (int128_t)q * f0 + (int128_t)r * g0;
    /* Verify that the bottom 62 bits of the result are zero, and then throw them away. */
    VERIFY_CHECK(((int64_t)cf & M62) == 0); cf >>= 62;
    VERIFY_CHECK(((int64_t)cg & M62) == 0); cg >>= 62;
    /* Compute limb 1 of t*[f,g], and store it as output limb 0 (= down shift). */
    cf += (int128_t)u * f1 + (int128_t)v * g1;
    cg += (int128_t)q * f1 + (int128_t)r * g1;
    f->v[0] = (int64_t)cf & M62; cf >>= 62;
    g->v[0] = (int64_t)cg & M62; cg >>= 62;
    /* Compute limb 2 of t*[f,g], and store it as output limb 1. */
    cf += (int128_t)u * f2 + (int128_t)v * g2;
    cg += (int128_t)q * f2 + (int128_t)r * g2;
    f->v[1] = (int64_t)cf & M62; cf >>= 62;
    g->v[1] = (int64_t)cg & M62; cg >>= 62;
    /* Compute limb 3 of t*[f,g], and store it as output limb 2. */
    cf += (int128_t)u * f3 + (int128_t)v * g3;
    cg += (int128_t)q * f3 + (int128_t)r * g3;
    f->v[2] = (int64_t)cf & M62; cf >>= 62;
    g->v[2] = (int64_t)cg & M62; cg >>= 62;
    /* Compute limb 4 of t*[f,g], and store it as output limb 3. */
    cf += (int128_t)u * f4 + (int128_t)v * g4;
    cg += (int128_t)q * f4 + (int128_t)r * g4;
    f->v[3] = (int64_t)cf & M62; cf >>= 62;
    g->v[3] = (int64_t)cg & M62; cg >>= 62;
    /* What remains is limb 5 of t*[f,g]; store it as output limb 4. */
    f->v[4] = (int64_t)cf;
    g->v[4] = (int64_t)cg;
}

/* Compute the inverse of x modulo modinfo->modulus, and replace x with it (constant time in x). */
static void secp256k1_modinv64(secp256k1_modinv64_signed62 *x, const secp256k1_modinv64_modinfo *modinfo) {
    /* Start with d=0, e=1, f=modulus, g=x, eta=-1. */
    secp256k1_modinv64_signed62 d = {{0, 0, 0, 0, 0}};
    secp256k1_modinv64_signed62 e = {{1, 0, 0, 0, 0}};
    secp256k1_modinv64_signed62 f = modinfo->modulus;
    secp256k1_modinv64_signed62 g = *x;
    int i;
    uint64_t eta = -(uint64_t)1;

    /* Do 12 iterations of 62 divsteps each = 744 divsteps. 724 suffices for 256-bit inputs. */
    for (i = 0; i < 12; ++i) {
        /* Compute transition matrix and new eta after 62 divsteps. */
        secp256k1_modinv64_trans2x2 t;
        eta = secp256k1_modinv64_divsteps_62(eta, f.v[0], g.v[0], &t);
        /* Update d,e using that transition matrix. */
        secp256k1_modinv64_update_de_62(&d, &e, &t, modinfo);
        /* Update f,g using that transition matrix. */
#ifdef VERIFY
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, -1) > 0); /* f > -modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, 1) <= 0); /* f <= modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &modinfo->modulus, -1) > 0); /* g > -modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &modinfo->modulus, 1) < 0);  /* g <  modulus */
#endif
        secp256k1_modinv64_update_fg_62(&f, &g, &t);
#ifdef VERIFY
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, -1) > 0); /* f > -modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, 1) <= 0); /* f <= modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &modinfo->modulus, -1) > 0); /* g > -modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &modinfo->modulus, 1) < 0);  /* g <  modulus */
#endif
    }

    /* At this point sufficient iterations have been performed that g must have reached 0
     * and (if g was not originally 0) f must now equal +/- GCD of the initial f, g
     * values i.e. +/- 1, and d now contains +/- the modular inverse. */
#ifdef VERIFY
    /* g == 0 */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &SECP256K1_SIGNED62_ONE, 0) == 0);
    /* |f| == 1, or (x == 0 and d == 0 and |f|=modulus) */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &SECP256K1_SIGNED62_ONE, -1) == 0 ||
                 secp256k1_modinv64_mul_cmp_62(&f, &SECP256K1_SIGNED62_ONE, 1) == 0 ||
                 (secp256k1_modinv64_mul_cmp_62(x, &SECP256K1_SIGNED62_ONE, 0) == 0 &&
                  secp256k1_modinv64_mul_cmp_62(&d, &SECP256K1_SIGNED62_ONE, 0) == 0 &&
                  (secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, 1) == 0 ||
                   secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, -1) == 0)));
#endif

    /* Optionally negate d, normalize to [0,modulus), and return it. */
    secp256k1_modinv64_normalize_62(&d, f.v[4] >> 63, modinfo);
    *x = d;
}

/* Compute the inverse of x modulo modinfo->modulus, and replace x with it (variable time). */
static void secp256k1_modinv64_var(secp256k1_modinv64_signed62 *x, const secp256k1_modinv64_modinfo *modinfo) {
    /* Start with d=0, e=1, f=modulus, g=x, eta=-1. */
    secp256k1_modinv64_signed62 d = {{0, 0, 0, 0, 0}};
    secp256k1_modinv64_signed62 e = {{1, 0, 0, 0, 0}};
    secp256k1_modinv64_signed62 f = modinfo->modulus;
    secp256k1_modinv64_signed62 g = *x;
    int i, j;
    uint64_t eta = -(uint64_t)1;
    int64_t cond;

    /* Do up to 12 iterations of 62 divsteps each = 744 divsteps, or until g=0 (whichever comes first). */
    for (i = 0; i < 12; ++i) {
        /* Compute transition matrix and new eta after 62 divsteps. */
        secp256k1_modinv64_trans2x2 t;
        eta = secp256k1_modinv64_divsteps_62_var(eta, f.v[0], g.v[0], &t);
        /* Update d,e using that transition matrix. */
        secp256k1_modinv64_update_de_62(&d, &e, &t, modinfo);
        /* Update f,g using that transition matrix. */
#ifdef VERIFY
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, -1) > 0); /* f > -modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, 1) <= 0); /* f <= modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &modinfo->modulus, -1) > 0); /* g > -modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &modinfo->modulus, 1) < 0);  /* g <  modulus */
#endif
        secp256k1_modinv64_update_fg_62(&f, &g, &t);
        /* If the bottom limb of g is zero, there is a chance that g=0. */
        if (g.v[0] == 0) {
            cond = 0;
            /* Check if the other limbs are also 0. */
            for (j = 1; j < 5; ++j) {
                cond |= g.v[j];
            }
            /* If so, we're done. */
            if (cond == 0) break;
        }
#ifdef VERIFY
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, -1) > 0); /* f > -modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, 1) <= 0); /* f <= modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &modinfo->modulus, -1) > 0); /* g > -modulus */
        VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &modinfo->modulus, 1) < 0);  /* g <  modulus */
#endif
    }

    /* At this point g is 0 and (if g was not originally 0) f must now equal +/- GCD of
     * the initial f, g values i.e. +/- 1, and d now contains +/- the modular inverse. */
#ifdef VERIFY
    /* g == 0 */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&g, &SECP256K1_SIGNED62_ONE, 0) == 0);
    /* |f| == 1, or (x == 0 and d == 0 and |f|=modulus) */
    VERIFY_CHECK(secp256k1_modinv64_mul_cmp_62(&f, &SECP256K1_SIGNED62_ONE, -1) == 0 ||
                 secp256k1_modinv64_mul_cmp_62(&f, &SECP256K1_SIGNED62_ONE, 1) == 0 ||
                 (secp256k1_modinv64_mul_cmp_62(x, &SECP256K1_SIGNED62_ONE, 0) == 0 &&
                  secp256k1_modinv64_mul_cmp_62(&d, &SECP256K1_SIGNED62_ONE, 0) == 0 &&
                  (secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, 1) == 0 ||
                   secp256k1_modinv64_mul_cmp_62(&f, &modinfo->modulus, -1) == 0)));
#endif

    /* Optionally negate d, normalize to [0,modulus), and return it. */
    secp256k1_modinv64_normalize_62(&d, f.v[4] >> 63, modinfo);
    *x = d;
}

#endif /* SECP256K1_MODINV64_IMPL_H */
