# GFMUL Implementation
These insights and implementation are mostly based on these three documents:
- **Intel Doc A**: (https://www.intel.com/content/dam/develop/external/us/en/documents/clmul-wp-rev-2-02-2014-04-20.pdf)
- **Intel Doc B**: (https://builders.intel.com/docs/networkbuilders/advanced-encryption-standard-galois-counter-mode-optimized-ghash-function-technology-guide-1693300747.pdf)
- **NIST**: (https://csrc.nist.rip/groups/ST/toolkit/BCM/documents/proposedmodes/gcm/gcm-revised-spec.pdf)
----
This document is structured into the following key points: 

0. [Mathematical notations](#mathematical-notations)
1. [Standard GMFUL implementation](#standard-gfmul-implementation)
2. [K-Optimization](#k-optimization)
3. [GHASH capable GFMUL implementation](#ghash-capable-gfmul-implementation)

----

## Mathematical Notations
GFMUL stands for Galois-Field multiplication. 
### GF($2^{128}$)
For AES-GCM this refers to the GF($2^{128}$) field.
Elements of this field can be expressed as polynoms with coefficients in GF($2^1$):
```math
\begin{align*}
a(x) \in \text{GF}(2^{128}) \
a(x) = \sum_{i=0}^{127}a_i x^i, \text{with } a_i \in \text{GF}(2^1)
\end{align*}
```
If the question now already arrives, how a `__m128i` value translates into the coefficients $a_i$?
--> E.g. whether:
```
__m128i a = _mm_set_epi32(0, 0, 0, 0x00000001);
```
results in:
```math
a(x) = 0x^{127} + 0x^{126} + ... + 0x^1 + 1x^0
```
or:
```math
a(x) = 1x^{127} + 0x^{126} + ... + 0x^0
```
This exactly boils down to the difference between the "classical/naive" GFMUL implementation and the GHASH necessary implementation.
Therefore, in order to stay consistent, we define that the polynom $b = 0x^{127}+...+0x^1+1x^0$ can be expressed with this variable 
```
__m128i b = _mm_set_epi32(0, 0, 0, 0x00000001);
```
### Addition and Multiplication
Addition and multiplication in GF(2) can be implemented as the **XOR** operation for addition, and the **AND** operation for multiplication.

Addition in GF(2^{128}) can be implemented as an bitwise **XOR** operation:
```math
\begin{align*}
a(x)+b(x)&=(a_{127}x^{127}+a_{126}x^{126}+...+a_0x^0) + (b_{127}x^{127}+b_{126}x^{126}+...+b_0x^0) \
&= (a_{127} \oplus b_{127})x^{127} + (a_{126}\oplus b_{126})x^{126} + ... + (a_0\oplus b_0)x^0 \
&= c(x)
\end{align*}
```
Multiplication can be implemented using first carry-less multiplication and then a reduction modulo the irreducible polynomial $P=x^{128}+x^7+x^2+x^1+1$.\
Multiplying two polynoms with degree $\leq 127$, leads to a new polynom with degree $\leq 254$.
This is not part of GF($2^{128}$) and must be reduced hence.

### Carry-Less Multiplication
This is not really mathematical background, more a combination of maths with computer science.\
Let's represent two polynomials $a(x)$ and $b(x)$ of degree $\leq 127$ using two polynomials each with degree $\leq 63$.
```math
\begin{align*}
a(x) = a_{1}\cdot x^{64} + a_{0} \
b(x) = b_{1}\cdot x^{64} + b_{0}
\end{align*}
```
This is because, we can not directly compute `CLMUL(a(x), b(x))`, but instead we can only multiply polynomials with degree $\leq 63$ (`$\text{CLMUL}_{64}$` --> `_mm_clmulepi64_si128()`)  yielding a polynomial of degree $\leq 255$.\
Thus:
```math
\text{CLMUL}(a(x)\text{, }b(x)) = (\text{CLMUL}_{64}(a_{1}(x), b_{1}(x)))\cdot x^{128} + (\text{CLMUL}_{64}(a_{0}(x), b_{1}(x)) + \text{CLMUL}_{64}(a_{1}(x), b_{0}(x)))\cdot x^{64} + (\text{CLMUL}_{64}(a_{0}(x), b_{0}(x)))
```

### Bit-Reflection
Because it is needed later on:
The bit-reflection transformation R is defined as follows for polynoms with degree $\leq 127$:
```math
R_{128}(a(x)) = x^{127}a(x^{-1})
```
Or in general:
```math
R_{n}(a(x)) = x^{n-1}a(x^{-1}), \text{if degree of a(x)} = n-1
```
If n is clear from context, we simply write $R$ or $(...)'$, e.g. $R(a(x)) = a(x)'$ .\
To give an example:
```math
\begin{align*}
a(x) &= x^{125} + x^{64} + x^{63} + x^1 + 1 \

R_{128}(a(x)) &= x^{127}\cdot a(x^{-1}) \
&= x^{127} \cdot (x^{-125} + x^{-64} + x^{-63} + x^{-1} + x^{-0}) \
&= x^{127-125} + x^{127-64} + x^{127-63} + x^{127-1} + x^{127-0} \
&= x^{2} + x^{63} + x^{64} + x^{126} + x^{127}
\end{align*}
```
In CS variable notation:
```
__m128i a      = _mm_set_epi32(0x20000000, 0x00000001, 0x80000000, 0x00000003);
__m128i a_refl = _mm_set_epi32(0xc0000000, 0x00000001, 0x80000000, 0x00000004);
```
Let's apply this on the irreducible polynomial P (will be useful later):
```math
\begin{align*}
P(x) &= x^{128} + x^7 + x^2 + x^1 + 1 \
R_{129}(Q(x)) &= x^{128} \cdot (x^{-128}+x^{-7}+x^{-2}+x^{-1}+x^{-0}) \
&= x^0+x^{121} + x^{126} + x^{127} + x^{128} \
&= x^{128} + x^{127}+x^{126}+x^{121}+1
\end{align*}
```

### CLMUL Identity
This following identity is crucial for implementing a correct gfmul implementation usable for GHASH.
```math
\text{CLMUL}(R_{128}(a(x))\text{, }R_{128}(b(x))) = R_{256}(\text{CLMUL}(a(x)\text{, }b(x)) << 1) 
```
Or in a more clean form:
```math
\text{CLMUL}(a(x)', b(x)') = (\text{CLMUL}(a(x), b(x)) << 1)'
```

## Standard GFMUL implementation
The standard (naive) gfmul implementation uses polynomials where the coefficient $a_0$ is stored at the least signifcant bit (LSB). \
Say $a(x)$ and $b(x)$, each polynomials inside GF($2^{128}$) with degree $\leq 127$, will be represented as two polynomials with degree $\leq 63} (64bit values):
```math
\begin{align*}
a(x) = a_{1}\cdot x^{64} + a_{0} \
b(x) = b_{1}\cdot x^{64} + b_{0}
\end{align*}
```
, and we have our irreducible polynomial:
```math
P(x) = x^{128} + x^7 + x^2 + x^1 + 1
```
Then:
```
GFMUL(a(x), b(x)) = CLMUL(a(x), b(x)) mod P(x) = c(x)
```
```math
\text{CLMUL}(a(x), b(x)) = (\text{CLMUL}_{64}(a_{1}(x), b_{1}(x)))\cdot x^{128} + (\text{CLMUL}_{64}(a_{0}(x), b_{1}(x)) + \text{CLMUL}_{64}(a_{1}(x), b_{0}(x)))\cdot x^{64} + (\text{CLMUL}_{64}(a_{0}(x), b_{0}(x))) = c(x)
```
or written more cleanly:
```math
\text{CLMUL}(a(x), b(x)) = (a_{1}b_{1}\cdot x^{128} + (a_{0}b_{1} + a_{1}b_{0})\cdot x^{64} + (a_{0}b_{0}) = c(x)
```
This results in a polynomial $c(x)$ which internally can be viewed as four smaller polynomials with degree $\leq 63$ each.
```math
c(x) = c_{3}x^{192} + c_{2}x^{128} + c_{1}x^{64} + c_{0}
```

If the result of the carry-less multiplication results in a polynomial with degree $\leq 127$ (Then $c_{3}$ and $c_{2}$ would each be 0 andn c(x) can be fitted in a single 128bit value), then no reduction is needed. \
Otherwise, the result can always be expressed partly as some form of:
```math
c(x) = x^{128}\dot (x^k + ...)
```
Because we know:
```math
x^{128} \equiv x^7 + x^2 + x^1 + 1 \mod P(x)
```
, we can replace the $x^{128}$ term with $(x^7+x^2+x^1+1)$ in $c(x)$, when calculating $c(x) \mod P(x)$.
We continue doing so, until its degree is $\leq 127$ such that $c(x)$ fits in GF($2^{128}$).

Intel Doc. B shows that this can be done by first calculating:
```math
u = c_{3} \cdot Q(x) = \text{CLMUL}_{64}(c_{3}, Q(x))
```
, with $Q(x) = x^7 + x^2 + x^1 + 1$.
Then:
```math
c*(x) = (c_{2}x^{128} + c_{1}x^{64} + c_{0}) + (u(x)x^{64})
```

## K-Optimization



# First Level Heading

Paragraph.

## Second Level Heading

Paragraph.

- bullet
+ other bullet
* another bullet
    * child bullet

1. ordered
2. next ordered

### Third Level Heading

Some *italic* and **bold** text and `inline code`.

An empty line starts a new paragraph.

Use two spaces at the end  
to force a line break.

A horizontal ruler follows:

---

Add links inline like [this link to the Qt homepage](https://www.qt.io),
or with a reference like [this other link to the Qt homepage][1].

    Add code blocks with
    four spaces at the front.

> A blockquote
> starts with >
>
> and has the same paragraph rules as normal text.

First Level Heading in Alternate Style
======================================

Paragraph.

Second Level Heading in Alternate Style
---------------------------------------

Paragraph.

[1]: https://www.qt.io

## GHASH capable GFMUL implementation
