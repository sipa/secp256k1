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
then evaluating $F(u, t)$. Encoding a given $x$ coordinate is done as follows:

* Loop:
  * Pick a random field element $u$.
  * Compute the set $L_u = F_{0,u}^{-1}(x)$ of $t$ values for which $F_{0,u}(t) = x$, which may have *0*, *4*, or *8* elements.
  * With probability $1 - \frac{\\#L_u}{8}$, restart the loop.
  * Select a random $t \in L_u$ and return $(u, t)$.

This is the *ElligatorSwift* algorithm, here given for just X coordinates. A version
that can encode full points will be given later.

## 2. The forward mapping function

First some definitions:
* Let the curve $E$ be the set of points $(x, y)$ for which $y^2 = x^3 + ax + b$, with $a$ and $b$
  public constants. For `secp256k1`, $a=0$ and $b=7$.
* Let the function $g(x) = x^3 + ax + b$, so the curve equation is also $y^2 = g(x)$.
* Define $V$ as the set of solutions $(x_1, x_2, x_3, z)$ to $z^2 = g(x_1)g(x_2)g(x_3)$.
* Define $S_u$ as the set of solutions $(X, Y)$ to $X^2 + (3u^2 + 4a)Y^2 = -g(u)$.

Note that for $V$, the the left hand side of the equation $z^2$ is a square, and thus the right
hand must also be square. Thus, out of the three right-hand side factors an even number (*0* or
*2*) must be non-squares. This implies that exactly *1* or exactly *3* out of
$\\{g(x_1), g(x_2), g(x_3)\\}$ must be a squares (unless $z=0$, but that's not possible for
curves of odd order like `secp256k1`). That means that for any valid $(x_1,x_2,x_3,z)$
on $V$, at least one of $\\{x_1, x_2, x_3\\}$ must be a valid X coordinate on $E$.

$F_{0,u}(t)$ will now be defined as:
* Let $(X, Y) = m_u(t)$, where $m_u$ is a function from field elements to $S_u$.
* Let $(x_1, x_2, x_3, z) = \psi_u(X, Y)$, where $\psi_u$ is a function from $S_u$ to $V$.
* Return the first element $x$ of $\[x_3, x_2, x_1\]$ (note the order!0 which is a valid X coordinate on $E$ (i.e., $g(x)$ is square).

The function $m_u$ depends on various properties of the curve, but is particularly simple for $a=0$
curves like `secp256k1`: $$m_u(t) = \left(\frac{u^3 + b + t^2}{2t}, \frac{u^3 + b - t^2}{2tu\sqrt{-3}}\right)$$
For $a \neq 0$ curves, see Equation 14 and Appendix A in the paper for $m_u$ (there defined as functions
$X(u,t)$ and $Y(u,t)$).

The function $\psi_u$ is the same for all curves: $\psi_u(X, Y) = (x_1, x_2, x_3, z)$ where:
* $x_1 = \frac{X}{2Y} - \frac{u}{2}$
* $x_2 = -\frac{X}{2Y} - \frac{u}{2}$
* $x_3 = u + 4Y^2$
* $z = \frac{g(x_3)(u^2 + ux_1 + x_1^2 + a)}{2Y}$

$$
\psi_u(X, Y) = (x_1, x_2, x_3, z) \text{where} \begin{itemize}
\item x_1 = \frac{X}{2Y} - \frac{u}{2}
\item x_2 = -\frac{X}{2Y} - \frac{u}{2}
\item x_3 = u + 4Y^2 \\ z = \frac{g(x_3)(u^2 + ux_1 + x_1^2 + a)}{2Y}
\end{itemize}
$$

