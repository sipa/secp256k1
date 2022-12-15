# Elligator Swift for secp256k1 explained

In this document we explain how the `ellswift` module implementation is related to the
construction in the
["SwiftEC: Shallue–van de Woestijne Indifferentiable Function To Elliptic Curves"](https://eprint.iacr.org/2022/759)
paper by Chavez-Saab, Rodrıguez-Henrıquez, and Tibouchi.

## 1. Introduction

The `ellswift` module effectively introduces a new 64-byte public key format, with the property
that (uniformly random) public keys can be encoded as 64-byte arrays which are computationally
indistinguishable from uniform byte arrays. The module provides functions to convert public keys
from and to this format, as well as convenience functions for key generation and ECDH that operate
directly on ellswift-encoded keys.

The encoding consists of the concatenation of two (32-byte big endian) encoded field elements $u$
and $t$. Together they encode an X coordinate on the curve $x$, or (see further) a full point on
the curve $Q = (x, y)$.

The mapping from $(u, t)$ to the curve is defined by the function $F(u, t) = F_{0,u}(t)$, where
$F_0$ will be defined later. Decoding is simply decoding the field elements $u$ and $t$, and
then evaluating $F(u, t)$. Encoding a given $x$ coordinate is conceptually done as follows:

* Loop:
  * Pick a uniformy random field element $u \neq 0$.
  * Compute the set $L = F_{0,u}^{-1}(x)$ of $t$ values for which $F_{0,u}(t) = x$, which may have *0*, *4*, or *8* elements.
  * With probability $1 - \frac{\\#L}{8}$, restart the loop.
  * Select a uniformly random $t \in L$ and return $(u, t)$.

This is the *ElligatorSwift* algorithm, here given for just X coordinates. An extension to full
$(x, y)$ points will be given later.

## 2. The decoding function

First some definitions:
* Let the curve $E$ be the set of points $(x, y)$ for which $y^2 = x^3 + ax + b$, with $a$ and $b$
  public constants.
  * For `secp256k1`, $a=0$ and $b=7$.
* Let the function $g(x) = x^3 + ax + b$, so the $E$ curve equation is also $y^2 = g(x)$.
* Define $V$ as the set of solutions $(x_1, x_2, x_3, z)$ to $z^2 = g(x_1)g(x_2)g(x_3)$.
* Define $S_u$ as the set of solutions $(X, Y)$ to $X^2 + (3u^2 + 4a)Y^2 = -g(u)$.
* $P_u$ is a function from field elements to $S_u$ that will be defined later.
* $\psi_u$ is a function from $S_u$ to $V$ that will be defined later.

Note that for $V$, the left hand side of the equation $z^2$ is a square, and thus the right
hand must also be square. Thus, out of the three right-hand side factors an even number (*0* or
*2*) must be non-squares. This implies that exactly *1* or exactly *3* out of
$\\{g(x_1), g(x_2), g(x_3)\\}$ must be square. That means that for any $(x_1,x_2,x_3,z) \in V$,
at least one of $\\{x_1, x_2, x_3\\}$ must be a valid X coordinate on $E$.

One exception to the above rule is the case where $z=0$, implying $g(x)=0$ for one of the $x_{\\{1,2,3\\}}$.
However, $g(x)=0$ implies $E$ has a point with Y coordinate $0$, which would have order 2. This is not
possible for any curve of odd order, like `secp256k1`.

$F_{0,u}(t)$ will now be defined as:
* Let $(X, Y) = P_u(t)$.
* Let $(x_1, x_2, x_3, z) = \psi_u(X, Y)$.
* Return the first element $x$ of $(x_3, x_2, x_1)$ (note the order!) which is a valid X coordinate on $E$ (i.e., $g(x)$ is square).

The function $P_u$ depends on various properties of the curve, but can be stated concisely for $a=0$
curves like `secp256k1`: $$P_u(t) = \left(\dfrac{u^3 + b - t^2}{2t}, \dfrac{u^3 + b + t^2}{2tu\sqrt{-3}}\right)$$
For $a \neq 0$ curves, see Equation 14 and Appendix A in the paper for $P_u(t)$ (there defined as functions
$X(u,t)$ and $Y(u,t)$ ).

The function $\psi_u$ is the same for all curves: $\psi_u(X, Y) = (x_1, x_2, x_3, z)$, where:

$$
\begin{array}{lcl}
  x_1 & = & \dfrac{X}{2Y} - \dfrac{u}{2} \\
  x_2 & = & -\dfrac{X}{2Y} - \dfrac{u}{2} \\
  x_3 & = & u + 4Y^2 \\
  z   & = & \dfrac{g(x_3)}{2Y}(u^2 + ux_1 + x_1^2 + a)
\end{array}
$$

### 2.1 Decoding for *a=0* curves

So put together, decoding $(u, t)$ to an X coordinate for $a=0$ curves corresponds to just:
* Let $X = \dfrac{u^3 + b - t^2}{2t}$.
* Let $Y = \dfrac{X + t}{u\sqrt{-3}}$.
* Return the first of $x$ in $(u + 4Y^2, \dfrac{-X}{2Y} - \dfrac{u}{2}, \dfrac{X}{2Y} - \dfrac{u}{2})$ for which $g(x)$ is square.

In the above procedure, there are 3 possible conditions that result in division by zero:
* $t=0$ yields a division by zero in the computation of $X$ and $Y$.
* $u=0$ yields a division by zero in the computation of $Y$.
* If $u^3 + b + t^2 = 0$, then $Y=0$, which results in a division by zero in the computation of $x_1$, $x_2$, and $z$.

In the paper, those divisions by zero all result in infinity being output due to the use of
projective coordinates. In the `ellswift` module, we avoid these special cases by making sure all
possible field elements $u$, $t$ map to valid X coordinates. To do so, the following approach is
used:
* Let $u'=u$ unless $u=0$, in which case $u'=1$ (guaranteeing $u' \neq 0$).
* Let $t'=t$ unless $t=0$, in which case $t'=1$ (guaranteeing $t' \neq 0$).
* Let $t''=t'$ unless $u'^3 + b + t'^2 = 0$, in which case $t''=2t'$ (guaranteeing $t'' \neq 0$ *and* $u'^3 + b + t''^2 \neq 0$).
* Run the normal decoding algorithm on $(u', t'')$:
  * Let $X = \dfrac{u'^3 + b - t''^2}{2t''}$.
  * Let $Y = \dfrac{X + t''}{u'\sqrt{-3}}$.
  * Return the first of $x$ in $(u' + 4Y^2, \dfrac{-X}{2Y} - \dfrac{u'}{2}, \dfrac{X}{2Y} - \dfrac{u'}{2})$ for which $g(x)$ is square.

This is implemented in `secp256k1_ellswift_xswiftec_frac_var` (which decodes to an X coordinate represented as a fraction).

## 3. The encoding function

To implement $F_{0,u}^{-1}(x)$, the function to find the set of inverses $t$ for which $F_{0,u}(t) = x$, we reverse the process:
* For each of the $x_{\\{1,2,3\\}}$ expressions in $\psi_u$, assume $x = x_i$ and solve for all $(X, Y)$ solutions on $S_u$.
  * For $x_1$ and $x_2$ there can be two such $(X, Y)$ solutions each.
  * For $x_3$ there can be four such $(X, Y)$ solutions.
* Compute the set of $t$ values for which $P_u(t) = (X, Y)$ for all these $(X, Y)$ solutions.
  * For $a=0$ curves, that set consists of $P_u^{-1}(X, Y) = u\sqrt{-3}Y - X$ applied to all $(X, Y)$.
* For each of the found $t$ values, verify that $F_{0,u}(t) = x$.
* Return the remaining $t$ values.

The third step here, verifying that $F_{0,u}(t) = x$, is necessary because for the $t$ values found through the $x_1$ and $x_2$ expressions,
it is possible that decoding through $x_3$ takes it to a valid X coordinate different from $x$,
which would take precedence over the $x_1$ or $x_2$ decoding. These $t$ values must be rejected.

Since we know that exactly one or exactly three out of $\\{x_1, x_2, x_3\\}$ are valid X coordinates for any $t$,
the case where either $x_1$ or $x_2$ is valid and in addition also $x_3$ is valid must mean that all three are valid.
This means that instead of checking whether $x_3$ is on the curve, it is also possible to check the other one out of
$x_1$ and $x_2$ is on the curve. This is significantly simpler, as it turns out.

Observe that $x_1 + x_2 = -u$. So given either $x = x_1$ or $x = x_2$, the other one of the two can be computed as
$-u - x$. Thus, when encoding $x$ through the $x_1$ or $x_2$ expressions, one can simply check whether $g(-u-x)$ is a square,
and if so, not include the corresponding $t$ values in the returned set. As this doesn't need $X$, $Y$, or $t$, this can in fact be computed
at the verify beginning of the algorithm.

It's not possible that an encoding found through the $x_1$ expression decodes to a different valid X coordinate using $x_2$ (which would
take precedence), for the same reason: if both $x_1$ and $x_2$ decodings were valid, $x_3$ would be valid as well, and thus take
precedence over both. Because of this, the $g(-u-x)$ being square test for $x_1$ and $x_2$ is the only test necessary to guarantee the found $t$
values round-trip back to the input $x$ correctly. This is the reason for choosing the $(x_3, x_2, x_1)$ precedence order in the decoder;
any other order requires more complicated round-trip checks in the encoder.

### 3.1 Switching to *v, w* coordinates

Before working out the formulas for all this, we switch to different variables for $S_u$. Let $v = \dfrac{X}{2Y} - \dfrac{u}{2}$, and
$w = 2Y$. Or in the other direction, $X = \dfrac{w}{2}(u+2v)$ and $Y = \dfrac{w}{2}$:
* $S_u'$ becomes the set of $(v, w)$ for which $w^2 (u^2 + uv + v^2 + a) = -g(u)$.
* We can restate $\psi_u$ in terms of $(v, w)$: $\psi_u'(v, w) = (x_1, x_2, x_3, z)$, where

$$
\begin{array}{lcl}
  x_1 & = & v \\
  x_2 & = & -u - v \\
  x_3 & = & u + w^2 \\
  z   & = & \dfrac{g(x_3)}{w}(u^2 + uv + v^2 + a)
\end{array}
$$

* $P_u^{-1}(X, Y) = u\sqrt{-3}Y - X$ can be stated for $(v,w)$ as $P_u^{'-1}(v, w) = w\left(\dfrac{\sqrt{-3}-1}{2}u - v\right)$.

We can now write the expressions for $\psi_u^{'-1}(x)$ explicitly, by solving each of the $x_{\\{1,2,3\\}}$
expressions for $v$ and $w$, also using the $S_u'$ equation $w^2 (u^2 + uv + v^2 + a) = -g(u)$ those
solutions have to satisfy:
* Assuming $x = x_1$, we find $v = x$ and $w = \pm\sqrt{-\dfrac{g(u)}{u^2 + uv + v^2 + a}}$.
* Assuming $x = x_2$, we find $v = -u-x$ and $w = \pm\sqrt{-\dfrac{g(u)}{u^2 + uv + v^2 + a}}$.
* Assuming $x = x_3$, we find $w = \pm\sqrt{x-u}$ and $v = \dfrac{-u}{2} \pm \dfrac{\sqrt{-w^2(4g(u) + w^2(3u^2 + 4a))}}{2w^2}$.

### 3.2 Avoiding computing all inverses

The ElligatorSwift algorithm as stated in Section 1 requires the computation of $L = F_{0,u}^{-1}(x)$ (the
set of all $t$ such that $(u, t)$ decode to $x$) in full. This is in fact not necessary.

Observe that the procedure of restarting with probability $1 - \dfrac{\\#L}{8}$ and otherwise returning a
uniformly random element from $L$ is actually equivalent to always padding $L$ with $\bot$ values up to length 8,
picking a uniformly random element from that, restarting whenever $\bot$ is picked:

* Loop:
  * Pick a uniformly random field element $u \neq 0$.
  * Compute the set $L = F_{0,u}^{-1}(x)$.
  * Let $T$ be the array consisting of the elements of $L$, plus $8 - \\#L$ times $\bot$.
  * Select a uniformly random $t \in T$.
  * If $t \neq \bot$, return $(u, t)$; restart loop otherwise.

Now notice that the order of elements in $T$ doesn't actually matter, as all we do is pick a uniformly
random element in it. As we have 8 distinct formulas for $\psi_u^{'-1}$ (taking the variants due to $\pm$ into account),
we can populate $T$ with the evaluations of those formulas, in order, taking into account:
* Formulas that yield no solutions, due to division by zero or taking the square root of a non-square, are made to return $\bot$.
* Any formula that returns $t=0$ is made to return $\bot$, as $t=0$ would be turned into $t'=1$ instead.
* For the $x_1$ and $x_2$ cases, if $g(-u-x)$ is a square, $\bot$ is returned instead.
* In case multiple formulas would return the same non- $\bot$ result, all but one of those must be turned into $\bot$ to avoid biasing those.

The last condition above only occurs with negligible probability for cryptographically-sized curves, but is interesting
to take into account as they allow exhaustive testing in small groups.

If we define $L = (C_{0,u}(x), C_{1,u}(x), C_{2,u}(x), \ldots, C_{7,u}(x))$, with each $C_{i,u}$ matching one of the formulas,
the loop can be simplified to only compute one of the inverses instead of all of them:

* Loop:
  * Pick a uniformly random field element $u \neq 0$.
  * Pick a uniformly random integer $c$ in $[0,8)$.
  * Let $t = C_{c,u}(x)$.
  * If $t \neq \bot$, return $(u, t)$; restart loop otherwise.

We map $c=0$ to the $x_1$ formula, $c=1$ to the $x_2$ formula, and $c=2$ and $c=3$ to the $x_3$ formula.
Those are then repeated as $c=4$ through $c=7$ for the other sign of $w$ (noting that in each formula, $w$ is a square root of some expression):

Define $C_{c,u}(x)$ as:
* If $c \\& 2$ (for $x_1$ and $x_2$ formulas):
  * If $g(-u-x)$ is square, return $\bot$ (as $x_3$ would be valid and take precedence).
  * If $c \\& 1$ (the $x_1$ formula) let $v = -u-x$, otherwise let $v = x$.
  * Let $s = \dfrac{-g(u)}{u^2 + uv + v^2 + a}$ (let $s = w^2$ in what follows).
* Otherwise (for $x_3$ formulas):
  * Let $s = x-u$.
  * Let $r = \sqrt{-s(4g(u) + s(3u^2+4a))};$ return $\bot$ if it doesn't exist.
  * Let $r' = -r$ if $c\,\\&\,1$, let $r' = r$ otherwise.
  * Let $v = \dfrac{r}{2s}-\dfrac{u}{2}$.
* Let $w = -\sqrt{s}$ if $c \\& 4$, $w = \sqrt{s}$ otherwise. Return $\bot$ if it doesn't exist.
* Return $P_u^{'-1}(v, w)$, which is $w\left(\dfrac{\sqrt{-3}-1}{2}u - v\right)$ for $a=0$ curves.
