# Elligator Swift for secp256k1 explained

In this document we explain how the `ellswift` module implementation is related to the
construction in the
["SwiftEC: Shallue–van de Woestijne Indifferentiable Function To Elliptic Curves"](https://eprint.iacr.org/2022/759)
paper by Chavez-Saab, Rodrıguez-Henrıquez, and Tibouchi.

* [1. Introduction](#1-introduction)
* [2. The decoding function](#2-the-decoding-function)
  + [2.1 Decoding for *a=0* curves](#21-decoding-for-a0-curves)
* [3. The encoding function](#3-the-encoding-function)
  + [3.1 Switching to *v, w* coordinates](#31-switching-to-v-w-coordinates)
  + [3.2 Avoiding computing all inverses](#32-avoiding-computing-all-inverses)
  + [3.3 Finding the inverse](#33-finding-the-inverse)
  + [3.4 Encoding for odd-ordered *a=0* curves](#34-encoding-for-odd-ordered-a0-curves)
* [4. Encoding and decoding full *(x, y)* coordinates](#4-encoding-and-decoding-full-x-y-coordinates)
  + [4.1 Full *(x, y)* coordinates for *a=0* curves](#41-full-x-y-coordinates-for-a0-curves)

## 1. Introduction

The `ellswift` module effectively introduces a new 64-byte public key format, with the property
that (uniformly random) public keys can be encoded as 64-byte arrays which are computationally
indistinguishable from uniform byte arrays. The module provides functions to convert public keys
from and to this format, as well as convenience functions for key generation and ECDH that operate
directly on ellswift-encoded keys.

The encoding consists of the concatenation of two (32-byte big endian) encoded field elements $u$
and $t$. Together they encode an x-coordinate on the curve $x$, or (see further) a full point $(x, y)$ on
the curve.

**Decoding** consists of decoding the field elements $u$ and $t$ (values above the field size $p$
are taken modulo $p$), and then evaluating $F_u(t)$, which for every $u$ and $t$ results in a valid
x-coordinate on the curve. The functions $F_u$ will be defined later.

**Encoding** a given $x$ coordinate is conceptually done as follows:
* Loop:
  * Pick a uniformy random field element $u$.
  * Compute the set $L = F_u^{-1}(x)$ of $t$ values for which $F_u(t) = x$, which may have up to *8* elements.
  * With probability $1 - \dfrac{\\#L}{8}$, restart the loop.
  * Select a uniformly random $t \in L$ and return $(u, t)$.

This is the *ElligatorSwift* algorithm, here given for just x-coordinates. An extension to full
$(x, y)$ points will be given later. The algorithm finds a uniformly random $(u, t)$ among (almost all) those
for which $F_u(t) = x$. Section 3.2 in the paper proves that the number of such encodings for
almost all x-coordinates on the curve (all but at most 39) is close to two times the field size
(specifically, it lies in the range $2q \pm (22\sqrt{q}-O(1))$, where $q$ is the size of the field).

## 2. The decoding function

First some definitions:
* $\mathbb{F}$ is the finite field of size $q$ (where $q \equiv 1 \Mod{3}$ and with characteristic 5 or higher).
  * For `secp256k1`, $q = 2^{256} - 2^{32} - 977$.
* Let the curve $E$ be the set of points $(x, y) \in \mathbb{F}^2$ for which $y^2 = x^3 + ax + b$, with $a$ and $b$
  public constants, that satisfy the conditions in Section 4.2 of the paper.
  * For `secp256k1`, $a=0$ and $b=7$.
* Let the function $g(x) = x^3 + ax + b$, so the $E$ curve equation is also $y^2 = g(x)$.
* Let the function $h(x) = 3x^3 + 4a$.
* Define $V$ as the set of solutions $(x_1, x_2, x_3, z)$ to $z^2 = g(x_1)g(x_2)g(x_3)$.
* Define $S_u$ as the set of solutions $(X, Y)$ to $X^2 + h(u)Y^2 = -g(u)$ and $Y \neq 0$.
* $P_u$ will be defined as function such that for almost all $(u, t)$ it holds that $P_u(t) \in S_u$.
* $\psi_u$ is a function from $S_u$ to $V$ that will be defined later.

**Note**: In the paper:
* $F_u$ corresponds to $F_{0,u}$ there.
* $P_u(t)$ corresponds to the pair of functions $X(u, t)$ and $Y(u, t)$ there, the result of which is jointly called $P$.
* All $S_u$ sets together correspond to $S$ there.
* All $\psi_u$ functions together (operating on elements of $S$) correspond to $\psi$ there.

Note that for $V$, the left hand side of the equation $z^2$ is square, and thus the right
hand must also be square. As the product of two non-squares is a square in $\mathbb{F}$,
out the three right-hand side factors an even number must be non-squares.
This implies that exactly *1* or exactly *3* out of
$\\{g(x_1), g(x_2), g(x_3)\\}$ must be square, and thus that for any $(x_1,x_2,x_3,z) \in V$,
at least one of $\\{x_1, x_2, x_3\\}$ must be a valid x-coordinate on $E$. There is one exception
to this, namely when $z=0$, but even then one of the three values is a valid x-coordinate.

**Define** the decoding function $F_u(t)$ as:
* Let $(x_1, x_2, x_3, z) = \psi_u(P_u(t))$.
* Return the first element $x$ of $(x_3, x_2, x_1)$ which is a valid x-coordinate on $E$ (i.e., $g(x)$ is square).

The function $P_u$ for $a \neq 0$ depends on values $X_0$ and $Y_0$, which are functions of $u$ and defined in Appendix A of the paper:

$$
P_u(t) = \left\\{\begin{array}{ll}
\left(\dfrac{g(u) - t^2}{2t} , \dfrac{g(u) + t^2}{2tu\sqrt{-3}}\right) & a = 0 \\
\left(\dfrac{g(u) + h(u)(Y_0 - tX_0)^2}{X_0(1+t^2h(u))} , \dfrac{tg(u) + (Y_0 - tX_0)(th(u)Y_0 + X_0)}{X_0(1+t^2h(u))}\right) & a \neq 0
\end{array}\right.
$$

It gives a valid result:
* For $a=0$, unless $u=0$, $t=0$, or $g(u) = -t^2$.
* For $a \neq 0$, unless $t^2h(u) = -1$ or $tg(u) = (tX_0 - Y_0)(th(u)Y_0 + X_0)$.

The function $\psi_u$ is the same for all curves: $\psi_u(X, Y) = (x_1, x_2, x_3, z)$, where:

$$
\begin{array}{lcl}
  x_1 & = & \dfrac{X}{2Y} - \dfrac{u}{2} && \\
  x_2 & = & -\dfrac{X}{2Y} - \dfrac{u}{2} && \\
  x_3 & = & u + 4Y^2 && \\
  z   & = & \dfrac{g(x_3)}{2Y}(u^2 + ux_1 + x_1^2 + a) \\
      & = & \dfrac{-g(u)g(x_3)}{8Y^3}
\end{array}
$$

### 2.1 Decoding for *a=0* curves

Put together, decoding $(u, t)$ to an x-coordinate for $a=0$ curves corresponds to just:

**Define** $F_u(t)$ as:
* Let $X = \dfrac{u^3 + b - t^2}{2t}$.
* Let $Y = \dfrac{X + t}{u\sqrt{-3}}$.
* Return the first $x$ in $(u + 4Y^2, \dfrac{-X}{2Y} - \dfrac{u}{2}, \dfrac{X}{2Y} - \dfrac{u}{2})$ for which $g(x)$ is square.

In the above procedure, there are 3 possible conditions that give undefined results:
* $t=0$ yields a division by zero in the computation of $X$ and $Y$.
* $u=0$ yields a division by zero in the computation of $Y$.
* If $u^3 + b + t^2 = 0$, then $Y=0$, which isn't valid on $S_u$.

**Note**: in the paper these conditions result in $\infty$ as output, due to the use of projective coordinates there.
We wish to avoid the need for callers to deal with this special case.

We make sure these conditions result in valid x-coordinates by altering the input in these cases:
* Let $u'=u$ unless $u=0$, in which case $u'=1$ (guaranteeing $u' \neq 0$).
* Let $t'=t$ unless $t=0$, in which case $t'=1$ (guaranteeing $t' \neq 0$).
* Let $t''=t'$ unless $u'^3 + b + t'^2 = 0$, in which case $t''=2t'$ (guaranteeing $t'' \neq 0$ *and* $u'^3 + b + t''^2 \neq 0$).
* Run the earlier algorithm on $(u', t'')$ instead of on $(u, t)$.

This is implemented in `secp256k1_ellswift_xswiftec_frac_var` (which decodes to an x-coordinate represented as a fraction), and
in `secp256k1_ellswift_xswiftec_var` (which outputs the actual x-coordinate).

## 3. The encoding function

To implement $F_u^{-1}(x)$, the function to find the set of inverses $t$ for which $F_u(t) = x$, we have to reverse the process:
* Find all the $(X, Y) \in S_u$ that could have given rise to $x$, through the $x_1$, $x_2$, or $x_3$ formulas in $\psi_u$.
* Map those $(X, Y)$ solutions to $t$ values using $P_u^{-1}(X, Y)$
  * For $a=0$ curves, $P_u^{-1}(X, Y) = u\sqrt{-3}Y - X$.
* For each of the found $t$ values, verify that $F_u(t) = x$.
* Return the remaining $t$ values.

The third step here, verifying that $F_u(t) = x$, is necessary because for the $(X, Y)$ values found through the $x_1$ and $x_2$ expressions,
it is possible that decoding through $\psi_u(X, Y)$ yields a valid $x_3$ on the curve, which would take precedence over the
$x_1$ or $x_2$ decoding. These $(X, Y)$ solutions must be rejected.

Since we know that exactly one or exactly three out of $\\{x_1, x_2, x_3\\}$ are valid x-coordinates for any $t$,
the case where either $x_1$ or $x_2$ is valid and in addition also $x_3$ is valid must mean that all three are valid.
This means that instead of checking whether $x_3$ is on the curve, it is also possible to check whether the other one out of
$x_1$ and $x_2$ is on the curve. This is significantly simpler, as it turns out.

Observe that $\psi_u$ guarantees that $x_1 + x_2 = -u$. So given either $x = x_1$ or $x = x_2$, the other one of the two can be computed as
$-u - x$. Thus, when encoding $x$ through the $x_1$ or $x_2$ expressions, one can simply check whether $g(-u-x)$ is a square,
and if so, not include the corresponding $t$ values in the returned set. As this doesn't need $X$, $Y$, or $t$, this can in fact be determined
before those values are computed.

It's not possible that an encoding found through the $x_1$ expression decodes to a different valid x-coordinate using $x_2$ (which would
take precedence), for the same reason: if both $x_1$ and $x_2$ decodings were valid, $x_3$ would be valid as well, and thus take
precedence over both. Because of this, the $g(-u-x)$ being square test for $x_1$ and $x_2$ is the only test necessary to guarantee the found $t$
values round-trip back to the input $x$ correctly. This is the reason for choosing the $(x_3, x_2, x_1)$ precedence order in the decoder;
any other order requires more complicated round-trip checks in the encoder.

### 3.1 Switching to *v, w* coordinates

Before working out the formulas for all this, we switch to different variables for $S_u$. Let $v = (X/Y - u)/2$, and
$w = 2Y$. Or in the other direction, $X = uw/2 + vw$ and $Y = w/2$:
* $S_u'$ becomes the set of $(v, w)$ for which $w^2 (u^2 + uv + v^2 + a) = -g(u)$ and $w \neq 0$.
* For $a=0$ curves, $P_u^{-1}$ can be stated for $(v,w)$ as $P_u^{'-1}(v, w) = w\left(\frac{\sqrt{-3}-1}{2}u - v\right)$.
* $\psi_u$ can be stated for $(v, w)$ as $\psi_u'(v, w) = (x_1, x_2, x_3, z)$, where

$$
\begin{array}{lcl}
  x_1 & = & v \\
  x_2 & = & -u - v \\
  x_3 & = & u + w^2 \\
  z   & = & \dfrac{g(x_3)}{w}(u^2 + uv + v^2 + a) \\
      & = & \dfrac{-g(u)g(x_3)}{w^3}
\end{array}
$$

We can now write the expressions for finding $(v, w)$ given $x$ explicitly, by solving each of the $\\{x_1, x_2, x_3\\}$
expressions for $v$ or $w$, and using the $S_u'$ equation to find the other variable:
* Assuming $x = x_1$, we find $v = x$ and $w = \pm\sqrt{-g(u)/(u^2 + uv + v^2 + a)}$.
* Assuming $x = x_2$, we find $v = -u-x$ and $w = \pm\sqrt{-g(u)/(u^2 + uv + v^2 + a)}$.
* Assuming $x = x_3$, we find $w = \pm\sqrt{x-u}$ and $v = -u/2 \pm \sqrt{-w^2(4g(u) + w^2h(u))}/(2w^2)$.

### 3.2 Avoiding computing all inverses

The ElligatorSwift algorithm as stated in Section 1 requires the computation of $L = F_u^{-1}(x)$ (the
set of all $t$ such that $(u, t)$ decode to $x$) in full. This is unnecessary.

Observe that the procedure of restarting with probability $(1 - \frac{\\#L}{8})$ and otherwise returning a
uniformly random element from $L$ is actually equivalent to always padding $L$ with $\bot$ values up to length 8,
picking a uniformly random element from that, restarting whenever $\bot$ is picked:

**Define** *ElligatorSwift(x)* as:
* Loop:
  * Pick a uniformly random field element $u$.
  * Compute the set $L = F_u^{-1}(x)$.
  * Let $T$ be the 8-element vector consisting of the elements of $L$, plus $8 - \\#L$ times $\\{\bot\\}$.
  * Select a uniformly random $t \in T$.
  * If $t \neq \bot$, return $(u, t)$; restart loop otherwise.

Now notice that the order of elements in $T$ doesn't matter, as all we do is pick a uniformly
random element in it, so we don't need to have all $\bot$ values at the end.
As we have 8 distinct formulas for finding $(v, w)$ (taking the variants due to $\pm$ into account),
we can associate every index in $T$ with exactly one of those formulas, making sure that:
* Formulas that yield no solutions, due to division by zero or taking the square root of a non-square, are made to return $\bot$.
* For the $x_1$ and $x_2$ cases, if $g(-u-x)$ is a square, $\bot$ is returned instead (the round-trip check).
* In case multiple formulas would return the same non- $\bot$ result, all but one of those must be turned into $\bot$ to avoid biasing those.

The last condition above only occurs with negligible probability for cryptographically-sized curves, but is interesting
to take into account as it allows exhaustive testing in small groups.

If we define $T = (G_{0,u}(x), G_{1,u}(x), \ldots, G_{7,u}(x))$, with each $G_{i,u}$ matching one of the formulas,
the loop can be simplified to only compute one of the inverses instead of all of them:

**Define** *ElligatorSwift(x)* as:
* Loop:
  * Pick a uniformly random field element $u$.
  * Pick a uniformly random integer $c$ in $[0,8)$.
  * Let $t = G_{c,u}(x)$.
  * If $t \neq \bot$, return $(u, t)$; restart loop otherwise.

This is implemented in `secp256k1_ellswift_xelligatorswift_var`.

### 3.3 Finding the inverse

To implement $G_{c,u}$, we map $c=0$ to the $x_1$ formula, $c=1$ to the $x_2$ formula, and $c=2$ and $c=3$ to the $x_3$ formula.
Those are then repeated as $c=4$ through $c=7$ for the other sign of $w$ (noting that in each formula, $w$ is a square root of some expression).
Ignoring the negligible cases, we get:

**Define** $G_{c,u}(x)$ as:
* If $c \in \\{0, 1, 4, 5\\}$ (for $x_1$ and $x_2$ formulas):
  * If $g(-u-x)$ is square, return $\bot$ (as $x_3$ would be valid and take precedence).
  * If $c \in \\{0, 4\\}$ (the $x_1$ formula) let $v = x$, otherwise let $v = -u-x$ (the $x_2$ formula)
  * Let $s = -g(u)/(u^2 + uv + v^2 + a)$ (using $s = w^2$ in what follows).
* Otherwise, when $c \in \\{2, 3, 6, 7\\}$ (for $x_3$ formulas):
  * Let $s = x-u$.
  * Let $r = \sqrt{-s(4g(u) + sh(u))}$.
  * Let $v = (r/s - u)/2$ if $c \in \\{3, 7\\}$; $(-r/s - u)/2$ otherwise.
* Let $w = \sqrt{s}$.
* Depending on $c$:
  * If $c \in \\{0, 1, 2, 3\\}$: return $P_u^{'-1}(v, w)$.
  * If $c \in \\{4, 5, 6, 7\\}$: return $P_u^{'-1}(v, -w)$.

Whenever a square root of a non-square is taken, $\bot$ is returned; for both square roots this happens with roughly
50% on random inputs. Similarly, when a division by 0 would occur, $\bot$ is returned as well; this will only happen
with negligible probability. On curves of odd order $u^2 + uv + v^2 + a \neq 0$, as equality would imply $g(-u-x) = g(x)$
which would mean the $g(-u-x)$ is square condition has triggered and $\bot$ would have been returned already
(however, on even-ordered curves it can also be a 0/0 as $g(u) = 0$ is possible).

**Note**: In the paper, the $case$ variable corresponds roughly to the $c$ above, but only takes on 4 possible values (1 to 4).
The conditional negation of $w$ at the end is done randomly, which is equivalent, but makes testing harder. We choose to
have the $G_{c,u}$ be deterministic, and capture all choices in $c$.

Now observe that the $c \in \\{1, 5\\}$ and $c \in \\{3, 7\\}$ conditions effectively perform the same $v \rightarrow -u-v$
transformation. Furthermore, that transformation has no effect on $s$ in the first branch
as $u^2 + ux + x^2 + a = u^2 + u(-u-x) + (-u-x)^2 + a$. Thus we can extract it out and move it down:

**Define** $G_{c,u}(x)$ as:
* If $c \in \\{0, 1, 4, 5\\}$:
  * If $g(-u-x)$ is square, return $\bot$.
  * Let $s = -g(u)/(u^2 + ux + x^2 + a)$.
  * Let $v = x$.
* Otherwise, when $c \in \\{2, 3, 6, 7\\}$:
  * Let $s = x-u$.
  * Let $v = (\sqrt{-s(4g(u) + sh(u))}/s - u)/2$.
* Let $w = \sqrt{s}$.
* Depending on $c$:
  * If $c \in \\{0, 2\\}$: return $P_u^{'-1}(v, w)$.
  * If $c \in \\{1, 3\\}$: return $P_u^{'-1}(-u-v, w)$.
  * If $c \in \\{4, 6\\}$: return $P_u^{'-1}(v, -w)$.
  * If $c \in \\{5, 7\\}$: return $P_u^{'-1}(-u-v, -w)$.

This shows there will always be exactly 0, 4, or 8 $t$ values for a given $(u, x)$ input.
There can be 0, 1, or 2 $(v, w)$ pairs before invoking $P_u^{'-1}$, and each results in 4 distinct $t$ values.

There are a few special conditions that only happen with negligible probability, and could optionally be dealt with:
* If $g(u) = 0$ then it is possible in the decoder that $g(x_1) = g(x_2) = 0$ while $g(x_3) \neq 0$, which violates the assumption
  that always either exactly 1 or exactly 3 out of $\\{x_1, x_2, x_3\\}$ are on the curve, which makes them hard to target by
  the encoder. Note that it's not necessary to outlaw $g(u)=0$ on the decoder side, as $g(x_2)$ will be on the curve.
  This condition cannot occur on curves of odd order.
* If $c \in \\{4,5,6,7\\}$ and $w=0$, return $\bot$ instead as otherwise the $c < 4$ and $c >= 4$ formulas would yield the same result.
  It doesn't need special care, as in the $c \in \\{4,5\\}$ case it's equivalent to $g(u)=0$ which is outlawed above, and
  in the $c \in \\{6,7\\}$ case it requires $s=0$ which is outlawed because of division by zero already.
* If $c \in \\{1,3,5,7\\}$ and $u + 2v = 0$, return $\bot$ instead as otherwise the odd and even $c$ expressions would return the same thing. This
  corresponds to either the two $x_3$ formulas yielding the same result, or the $x_1$ and $x_2$ formulas yielding the same result.
  The latter condition is impossible because it implies $g(-u-x) = g(x)$, which means the roundtrip test would have failed already.
* If $c \in \\{0, 1, 4, 5\\}$ and $g(x)=0$, return $\bot$ instead as otherwise the $x_3$ formulas would be equal to the $x_1$ or $x_2$ formulas,
  and we have to prefer the $x_3$ formulas. This condition is not possible for odd-ordered curves.
* Any input for which $P_u^{'-1}$ doesn't exist has to return $\bot$; for $a=0$ curves this happens when $u=0$.

### 3.4 Encoding for odd-ordered *a=0* curves

Specialized for odd-ordered $a=0$ curves:

**Define** $G_{c,u}(x)$ as:
* If $u=0$, return $\bot$.
* If $c \in \\{0, 1, 4, 5\\}$:
  * If $g(-u-x)$ is square, return $\bot$
  * Let $s = -g(u)/(u^2 + ux + x^2 + a)$ (note: can't cause division by 0).
  * Let $v = x$.
* Otherwise, when $c \in \\{2, 3, 6, 7\\}$:
  * Let $s = x-u$.
  * Let $v = (\sqrt{-s(4g(u) + sh(u))}/s - u)/2$.
  * If $c \in \\{3, 7\\}$ and $u + 2v = 0$, return $\bot$.
* Let $w = \sqrt{s}$.
* Depending on $c$:
  * If $c \in \\{0, 2\\}$: return $w(\frac{\sqrt{-3}-1}{2}u - v)$.
  * If $c \in \\{1, 3\\}$: return $w(\frac{\sqrt{-3}+1}{2}u + v)$.
  * If $c \in \\{4, 6\\}$: return $w(\frac{-\sqrt{-3}+1}{2}u + v)$.
  * If $c \in \\{5, 7\\}$: return $w(\frac{-\sqrt{-3}-1}{2}u - v)$.

This is implemented in `secp256k1_ellswift_xswiftec_inv_var`.

And the x-only ElligatorSwift encoding algorithm is still:

**Define** *ElligatorSwift(x)* as:
* Loop:
  * Pick a uniformly random field element $u$.
  * Pick a uniformly random integer $c$ in $[0,8)$.
  * Let $t = G_{c,u}(x)$.
  * If $t \neq \bot$, return $(u, t)$; restart loop otherwise.

For the targetted curves, the $T = (G_{0,u}(x), \ldots, G_{7,u}(x))$ induced by the above
$G_{c,u}$ implementation has the property that it contains every $t$ for which $F_u(t) = x$
*exactly* once, with the exception of the remapped $u=0$, $t=0$, and $u^3 + b + t^2 = 0$
cases, which happen with negligible probability.

## 4. Encoding and decoding full *(x, y)* coordinates

So far we have only addressed encoding and decoding x-coordinates, but in some cases an encoding
for full points with $(x, y)$ coordinates is desirable. It is possible to encode this information
in $t$ as well.

Note that for any $(X, Y) \in S_u$, $(\pm X, \pm Y)$ are all on $S_u$. Moreover, all of these are
mapped to the same x-coordinate. Negating $X$ or negating $Y$ just results in $x_1$ and $x_2$
being swapped, and doesn't affect $x_3$. This won't change the outcome x-coordinate as the order
of $x_1$ and $x_2$ only matters if both were to be valid, in which case $x_3$ would be used instead.

Still, these four $(X, Y)$ combinations all correspond to distinct $t$ values, so we can encode
the sign of the y-coordinate in the sign of $X$ or the sign of $Y$. They correspond to the
four distinct $P_u^{'-1}$ calls in the definition of $G_{u,c}$.

**Note**: In the paper, the sign of the y coordinate is encoded in a separately-coded bit.

To encode the sign of $y$ in the sign of $Y$:

**Define** *Decode(u, t)* for full $(x, y)$ as:
* Let $(X, Y) = P_u(t)$.
* Let $x$ be the first value in $(u + 4Y^2, \frac{-X}{2Y} - \frac{u}{2}, \frac{X}{2Y} - \frac{u}{2})$ for which $g(x)$ is square.
* Let $y = \sqrt{g(x)}$.
* If $sign(y) = sign(Y)$, return $(x, y)$; otherwise return $(x, -y)$.

And encoding would be done using a $G_{c,u}(x, y)$ function defined as:

**Define** $G_{c,u}(x, y)$ as:
* If $c \in \\{0, 1\\}$:
  * If $g(-u-x)$ is square, return $\bot$
  * Let $s = -g(u)/(u^2 + ux + x^2 + a)$.
  * Let $v = x$.
* Otherwise, when $c \in \\{2, 3\\}$:
  * Let $s = x-u$.
  * Let $v = (\sqrt{-s(4g(u) + sh(u))}/s - u)/2$.
* Let $w = \sqrt{s}$.
* Let $w' = w$ if $sign(w/2) = sign(y)$; otherwise $w' = -w$.
* Depending on $c$:
  * If $c \in \\{0, 2\\}$: return $P_u^{'-1}(v, w')$.
  * If $c \in \\{1, 3\\}$: return $P_u^{'-1}(-u-v, w')$.

Note that $c$ now only ranges $[0,4)$, as the sign of $w'$ is decided based on that of $y$, rather than on $c$.

In the above $sign$ could be implemented in several ways, such as parity of the integer representation
of the input field element (for prime-sized fields) or the quadratic residuosity (for fields where
$-1$ is not square). The choice doesn't matter, as long as $sign(0)=0$, and $1$ or $-1$ for every other
value, and has the property that $-sign(x) = sign(-x)$.

### 4.1 Full *(x, y)* coordinates for *a=0* curves

For $a=0$ curves, there is another option. Note that for those,
the $P_u(t)$ function translates negations of $t$ to negations of (both) $X$ and $Y$. Thus, we can use $sign(t)$ to
encode the y-coordinate. Combined with the earlier remapping to guarantee all inputs land on the curve, we get
as decoder:

**Define** *Decode(u, t)* as:
* Let $u'=u$ unless $u=0$, in which case $u'=1$.
* Let $t'=t$ unless $t=0$, in which case $t'=1$.
* Let $t''=t'$ unless $u'^3 + b + t'^2 = 0$, in which case $t''=2t'$.
* Let $X = \dfrac{u'^3 + b - t''^2}{2t''}$.
* Let $Y = \dfrac{X + t''}{u'\sqrt{-3}}$.
* Let $x$ be the first element of $(u' + 4Y^2, \frac{-X}{2Y} - \frac{u'}{2}, \frac{X}{2Y} - \frac{u'}{2})$ for which $g(x)$ is square.
* Let $y = \sqrt{g(x)}$.
* If $sign(y) = sign(t)$, return $(x, y)$; otherwise return $(x, -y)$.

This is implemented in `secp256k1_ellswift_swiftec_var`. The used $sign(x)$ function is the parity of $x$ when represented as in integer in $[0,q)$.

The corresponding encoder would invoke the x-only one, but negating the output $t$ if $sign(t) \neq sign(y)$.

This is implemented in `secp256k1_ellswift_elligatorswift_var`.
