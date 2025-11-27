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
a(x) \in \text{GF}(2^{128}) \implies a(x) := \sum_{i=0}^{127}a_i x^i 
\end{align*}
```
, with $a_i \in$ GF($2^1$) 

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
Therefore, in order to stay consistent for the "normal" approach, we define that the polynom $b = 0x^{127}+...+0x^1+1x^0$ can be expressed with this variable 
```
__m128i b = _mm_set_epi32(0, 0, 0, 0x00000001);
```
### Addition and Multiplication
Addition and multiplication in GF(2) can be implemented as the **XOR** operation for addition, and the **AND** operation for multiplication.

Addition in GF(2^{128}) can be implemented as an bitwise **XOR** operation:
```math
\begin{align*}
a(x)+b(x)&=\big(a_{127}x^{127}+a_{126}x^{126}+...+a_0x^0\big) + \big(b_{127}x^{127}+b_{126}x^{126}+...+b_0x^0\big) \
&= (a_{127} \oplus b_{127})x^{127} + (a_{126}\oplus b_{126})x^{126} + ... + (a_0\oplus b_0)x^0 \
&= c(x)
\end{align*}
```
Multiplication can be implemented using first carry-less multiplication and then a reduction with modulo the irreducible polynomial $P(x) := x^{128}+x^7+x^2+x^1+1$.\
Multiplying two polynoms with degree $\leq 127$, leads to a new polynom with degree $\leq 254$.
This is not part of GF($2^{128}$) and must be reduced hence.

### Splitting polynomials
This routine and technique is important for the following algorithm. \
Whenever we have a polynomial $f(x)$ of degree $n-1$, we can split it up into $p$ smaller equally sized polynomials $`f_{[j]}(x)`$, if $p|n$ (p is divisor of n), where each smaller polynomial is of max degree $k = \frac{n}{p}$.
```math
f(x) := \sum_{j=0}^{p}{f_{[j]}(x)x^{k\cdot j}}
```
Examples:
- For a polynomial $a$ of degree 127 ($n=128$), being a standard element of GF($2^{128}$):
```math
a(x) = a_{[1]}(x)x^{64} + a_{[0]}(x)
```
   where $a_{\[1\]}$ and $a_{[0]}$ are both polynomials of degree $\leq 63$.
- For a polynomial $b$ of degree 255 ($n=256$), we can split it up the following way:
```math
b(x) = b_{[3]}(x)x^{192} + b_{[2]}(x)x^{128} + b_{[1]}(x)x^{64} + b_{[0]}(x)
```
<br>
(If you think about it, this is nothing else than in order to store a polynom of degree 127 in memory, we store it in two 64bit registers and just keep track of which one is the higher and which one is the lower.)

### Carry-Less Multiplication
This is not really mathematical background, more a combination of maths with computer science.\
Let's represent two polynomials $a(x)$ and $b(x)$ of degree $\leq 127$ using two polynomials each with degree $\leq 63$.
```math
\begin{align*}
a(x) := a_{[1]}(x)x^{64} + a_{[0]}(x) \
b(x) := b_{[1]}(x)x^{64} + b_{[0]}(x)
\end{align*}
```
This is because, we can not directly compute `CLMUL(a(x), b(x))`, but instead we can only multiply polynomials with degree $\leq 63$ ($`\text{CLMUL}_{64}\rightarrow`$ `_mm_clmulepi64_si128()`)  yielding a polynomial of degree $\leq 255$.\
Thus:
```math
\begin{align*}
\text{CLMUL}(a(x), b(x)) &= \big(\text{CLMUL}_{64}(a_{[1]}(x), b_{[1]}(x))\big) x^{128} \
&+ \big(\text{CLMUL}_{64}(a_{[0]}(x), b_{[1]}(x)) + \text{CLMUL}_{64}(a_{[1]}(x), b_{[0]}(x))\big) x^{64} \
&+ \big(\text{CLMUL}_{64}(a_{[0]}(x), b_{[0]}(x))\big)
\end{align*}
```
or written more cleanly (but not mathematically correct):
```math
\text{CLMUL}(a(x), b(x)) = (a_{[1]}b_{[1]})\cdot x^{128} + (a_{[0]}b_{[1]} + a_{[1]}b_{[0]}) x^{64} + (a_{[0]}b_{[0]})
```

### Bit-Reflection
Because it is needed later on:
The bit-reflection transformation R is defined as follows for polynoms with degree $= 127$:
```math
R_{128}(a(x)) := x^{127}\cdot a(x^{-1})
```
Or in general:
```math
R_{n}(a(x)) := x^{n-1}\cdot a(x^{-1}), \text{if degree of a(x)} = n-1
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
Let's apply this on the irreducible polynomial P $(n=129)$ (will be useful later):
```math
\begin{align*}
P(x) &:= x^{128} + x^7 + x^2 + x^1 + 1 \\[6pt]
R_{129}(P(x)) &= x^{128} \cdot (x^{-128}+x^{-7}+x^{-2}+x^{-1}+x^{-0}) \
&= x^0+x^{121} + x^{126} + x^{127} + x^{128} \
&= x^{128} + x^{127}+x^{126}+x^{121}+1
\end{align*}
```

### CLMUL Identity
This following identity is crucial for implementing a correct gfmul implementation usable for GHASH.
```math
\text{CLMUL}\big(R_{128}(a(x)), R_{128}(b(x))\big) = R_{256}\big(\text{CLMUL}(a(x), b(x)) << 1\big) 
```
Or in a more clean form:
```math
\text{CLMUL}\big(a(x)', b(x)'\big) = \big(\text{CLMUL}(a(x), b(x)) << 1\big)'
```

## Standard GFMUL implementation
The standard (naive) gfmul implementation uses polynomials where the coefficient $a_0$ is stored at the least signifcant bit (LSB). \
Say $a(x)$ and $b(x)$, each polynomials inside GF($2^{128}$) with degree $\leq 127$, will be splitted into two polynomials with degree $\leq 63$ (64bit values):
```math
\begin{align*}
a(x) := a_{[1]}(x)x^{64} + a_{[0]}(x) \
b(x) := b_{[1]}(x)x^{64} + b_{[0]}(x)
\end{align*}
```
, and we have our irreducible polynomial and $Q(x)$:
```math
\begin{align*}
P(x) &:= x^{128} + x^7 + x^2 + x^1 + 1 \
Q(x) &:= x^7 + x^2 + x^1 + 1
\end{align*}
```
Then:
```
GFMUL(a(x), b(x)) = CLMUL(a(x), b(x)) mod P(x) = c(x)
```
```math
\begin{align*}
c(x) &:= \text{CLMUL}(a(x), b(x))  \\[4pt]
&= \big(\text{CLMUL}_{64}(a_{[1]}(x), b_{[1]}(x))\big)x^{128} \
&+ \big(\text{CLMUL}_{64}(a_{[0]}(x), b_{[1]}(x)) + \text{CLMUL}_{64}(a_{[1]}(x), b_{[0]}(x))\big)x^{64} \
&+ \big(\text{CLMUL}_{64}(a_{[0]}(x), b_{[0]}(x))\big)
\end{align*}
```
or written more cleanly:
```math
c(x) \approx (a_{[1]}b_{[1]})x^{128} + \big(a_{[0]}b_{[1]} + a_{[1]}b_{[0]}\big)x^{64} + (a_{[0]}b_{[0]}) = c(x)
```
This results in a polynomial $c(x)$ which internally can be viewed as four smaller polynomials with degree $\leq 63$ each (keep in mind that $a_{\[1\]}b_{\[1\]}$ is of degree $\leq 127$.
```math
\begin{align*}
c(x) &= c_{[3]}(x)x^{192} + c_{[2]}(x)x^{128} + c_{[1]}(x)x^{64} + c_{[0]}(x) \
&\approx c_{[3]}x^{192} + c_{[2]}x^{128} + c_{[1]}x^{64} + c_{[0]}
\end{align*}
```

If the result of the carry-less multiplication results in a polynomial with degree $\leq 127$ (Then $c_{[3]}$ and $c_{[2]}$ would each be 0 and c(x) can be fitted in a single 128bit value), then no reduction is needed. \
Otherwise, we need to calculate $c(x) \mod P(x)$:
We use:
```math
x^{128} \equiv x^7 + x^2 + x^1 + 1 \mod \big(x^{128} + x^7 + x^2 + x^1 + 1 = P(x)\big)
```
We do the reduction starting with the polynomial $c(x)$ of degree $\leq 254$:
```math
\begin{align*}
c(x) &\equiv c_{[3]}(x)x^{64}x^{128} + c_{[2]}(x)x^{128} + c_{[1]}(x)x^{64} + c_{[0]}(x) \mod P(x) \\[6pt]
&\equiv \big(c_{[3]}(x) \cdot (x^7+x^2+x^1+1)\big)x^{64} + c_{[2]}(x)x^{128} + c_{[1]}(x)x^{64} + c_{[0]}(x) \mod P(x) \\[6pt]
\end{align*}
```
We can write:
```math
\begin{align*}
u(x) &:= (c_{[3]}(x) \cdot (x^7+x^2+x^1+1)) \\[6pt]
&= \text{CLMUL}_{64}(c_{[3]}(x), Q(x)) \\[6pt]
&= u_{[1]}(x)x^{64} + u_{[0]}(x) \\[6pt]
\end{align*}
```
We split $u(x)$ again in two smaller polynomials. \
Hence:
```math
\begin{align*}
c(x) &\equiv \big(c_{[3]}(x) \cdot (x^7+x^2+x^1+1)\big)x^{64} + c_{[2]}(x)x^{128} + c_{[1]}(x)x^{64} + c_{[0]}(x) \mod P(x) \\[6pt]
&\equiv \big(u_{[1]}(x)x^{64} + u_{[0]}(x)\big)x^{64} + c_{[2]}(x)x^{128} + c_{[1]}(x)x^{64} + c_{[0]}(x) \mod P(x)  \\[6pt]
&\equiv \big(u_{[1]}(x) + c_{[2]}(x)\big)x^{128} + \big(u_{[0]}(x) + c_{[1]}(x)\big)x^{64} + c_{[0]}(x) \mod P \\[6pt]
&\equiv \big((u_{[1]}(x) + c_{[2]}(x)) \cdot (x^7+x^2+x^1+1)\big) + \big(u_{[0]}(x) + c_{[1]}(x)\big)x^{64} + c_{[0]}(x) \mod P \\[6pt]
\end{align*}
```
We substitute again:
```math
\begin{align*}
v(x) &:= (u_{[1]}(x) + c_{[2]}(x)) \cdot (x^7+x^2+x^1+1) \\[6pt]
&= \text{CLMUL}_{64}\big((u_{[1]}(x) + c_{[2]}(x)), Q(x)\big) \\[6pt]
&= v_{[1]}(x)x^{64} + v_{[0]}(x)  \\[6pt]
\end{align*}
```
We split $v(x)$ again in two smaller polynomials. \
Combined:
```math
\begin{align*}
c(x) &\equiv \big((u_{[1]}(x) + c_{[2]}(x)) \cdot (x^7+x^2+x^1+1)\big) + \big(u_{[0]}(x) + c_{[1]}(x)\big)x^{64} + c_{[0]}(x) \mod P \\[6pt]
&\equiv v_{[1]}(x)x^{64} + v_{[0]}(x) + \big(u_{[0]}(x) + c_{[1]}(x)\big)x^{64} + c_{[0]}(x) \mod P \\[6pt]
&\equiv \big(v_{[1]}(x) + u_{[0]}(x) + c_{[1]}(x)\big)x^{64} + (v_{[0]}(x) + c_{[0]}(x)) \\[6pt]
\end{align*}
```

The complete sketch can also be viewed here:
<img width="942" height="1234" alt="grafik" src="https://github.com/user-attachments/assets/4daf5663-f696-4a7b-8ddb-3b41ef948c61" />

Image is from Doc. B.
Please note difference in notations.

### Code Implementation
See a full code implementation below using Intel AVX instruction set
```c
__m128i gfmul(__m128i a, __m128i b){
    __m128i q = _mm_set_epi32(0, 0, 0, 0x00000087);   // this corresponds to Q(x) = x^7 + x^2 + x^1 + 1

    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);

    __m128i mid = _mm_xor_si128(a0b1, a1b0); // computes mid = A0B1 + A1B0

    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8)); // computes C[1:0] = A0B0 + (mid << x^64)
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8)); // computes C[3:2] = A1B1 + (mid >> x^64)

    // Step 2.1: Reduce
    __m128i x = _mm_clmulepi64_si128(a1b1, q, 0x01);  // computes C[3] * Q = upper(A1B1) * lower(Q)
    c01 = _mm_xor_si128(c01, _mm_slli_si128(x, 8)); // add lower half of x (X[0]) to upper part of C[1:0]
    c23 = _mm_xor_si128(c23, _mm_srli_si128(x, 8)); // add higher half of x (X[1]) to lower part of C[3:2] (higher part is just dangling around)

    // Step 2.2: Reduce
    x = _mm_clmulepi64_si128(c23, q, 0x00);         // works because higher part is not used in the calculation
    c01 = _mm_xor_si128(c01, x);

    return c01;
}
```

## K-Optimization
Intel Doc. B shows how we can optimize this routine by looking at the CLMUL calculation result just before reducing with $P(x)$:
```math
\begin{align*}
c(x) &\equiv \big(a_{[1]}(x)b_{[1]}(x)\big)x^{128} + \big(a_{[0]}(x)b_{[1]}(x) + a_{[1]}(x)b_{[0]}(x)\big)x^{64} + (a_{[0]}(x)b_{[0]}(x)) \mod P(x) \\[6pt]
&\equiv \big(a_{[1]}(x)\cdot b_{[1]}(x)\cdot Q(x)\big) + \big(a_{[0]}(x)b_{[1]}(x) + a_{[1]}(x)b_{[0]}(x)\big)x^{64} + (a_{[0]}(x)b_{[0]}(x)) \mod P(x) \\[6pt]
\end{align*}
```
The K-optimization is realized through knowing that in a AES-GCM execution, the encryption key does not change, and with that the GHASH key $H$ does not change.
The routine of calculating the GHASH requires multiplying a dynamic value ($a$) with fixed variations of $H$ (e.g. $H^1, H^2, H^3, ...$).
These values can be pre-computed and thus $b$ for multiplying is always known beforehand.

Let's define $K$:
```math
\begin{align*}
K(x) &:= \text{CLMUL}_{64}(b_{[1]}(x), Q(x)) \\
&= K_{[1]}(x)x^{64} + K_{[0]}(x)
\end{align*} 
```
then:
```math
\begin{align*}
c(x) &\equiv \big(a_{[1]}(x)\cdot b_{[1]}(x)\cdot Q(x)\big) + \big(a_{[0]}(x)b_{[1]}(x) + a_{[1]}(x)b_{[0]}(x)\big)x^{64} + (a_{[0]}(x)b_{[0]}(x)) \mod P(x) \\[6pt]
&\equiv \big(a_{[1]}(x)\cdot K(x)\big) + \big(a_{[0]}(x)b_{[1]}(x) + a_{[1]}(x)b_{[0]}(x)\big)x^{64} + (a_{[0]}(x)b_{[0]}(x)) \mod P(x) \\[6pt]
&\equiv \big(a_{[1]}(x)\cdot (K_{[1]}(x)x^{64} + K_{[0]}(x))\big) + \big(a_{[0]}(x)b_{[1]}(x) + a_{[1]}(x)b_{[0]}(x)\big)x^{64} + (a_{[0]}(x)b_{[0]}(x)) \mod P(x) \\[6pt]
&\equiv (a_{[1]}(x)K_{[1]}(x))x^{64} + (a_{[1]}(x)K_{[0]}(x)) + \big(a_{[0]}(x)b_{[1]}(x) + a_{[1]}(x)b_{[0]}(x)\big)x^{64} + (a_{[0]}(x)b_{[0]}(x)) \mod P(x) \\[6pt]
\end{align*}
```
With that, $c(x)$ has already been reduced such that the polynomial's degree $\leq 192$, and thus only one more reduction is required.
```math
c(x) \equiv c_{[2]}(x)x^{128} + c_{[1]}(x)x^{64} + c_{[0]}(x) \mod P(x)
```
Additionally, since:
```math
\begin{align*}
c(x) &\equiv \big((a_{[1]}(x)K_{[1]}(x))x^{64} + (a_{[1]}(x)K_{[0]}(x))\big) + \big(a_{[0]}(x)b_{[1]}(x) + a_{[1]}(x)b_{[0]}(x)\big)x^{64} + (a_{[0]}(x)b_{[0]}(x)) \mod P(x) \\[6pt]
&\equiv \big(a_{[1]}(x)K_{[1]}(x) + a_{[0]}(x)b_{[1]}(x) + a_{[1]}(x)b_{[0]}(x)\big)x^{64} + (a_{[1]}(x)K_{[0]}(x)) + (a_{[0]}(x)b_{[0]}(x)) \mod P(x) \\[6pt]
&\equiv \big(a_{[1]}(x)\cdot (K_{[1]}(x) + b_{[0]}(x))\big)x^{64} + (a_{[0]}(x)b_{[1]}(x))x^{64} + (a_{[1]}(x)K_{[0]}(x)) + (a_{[0]}(x)b_{[0]}(x)) \mod P(x) \\[6pt]
\end{align*}
```
Hence, if we update $K$ to:
```math
\begin{align*}
K^{*}(x) &:= \text{CLMUL}_{64}(b_{[1]}(x), Q(x)) + b_{[0]}(x)x^{64} \\[4pt]
&= K^{*}_{[1]}(x)x^{64} + K^{*}_{[0]}(x)
\end{align*}
```
the equation is further reduced to:
```math
c(x) \equiv \big((a_{[1]}(x)K^{*}_{[1]}(x)) + (a_{[0]}(x)b_{[1]}(x))\big)x^{64} + (a_{[1]}(x)K^{*}_{[0]}(x)) + (a_{[0]}(x)b_{[0]}(x)) \mod P(x)
```
While still requiring one reduction, it also requires one less $\text{CLMUL}_{64}$ invocation by extended pre-computing of $K^{*}(x)$.

This complete sketch can be viewed here:
<img width="1040" height="1208" alt="grafik" src="https://github.com/user-attachments/assets/16216fce-53ff-4335-b60b-056e7838bcbd" />

Again this image is also from Doc. B. 
Please note the different notation scheme.
The $K$ used in the image corresponds to our $K^{*}(x)$.

### Code Implementation
See a full code implementation below using Intel AVX instruction set which calculates $K^{*}(x)$ (`k`) every invocation.
```c
__m128i gfmul_k_optimized(__m128i a, __m128i b){
    // Step 0: Pre-requesites
    __m128i q = _mm_set_epi32(0, 0, 0, 0x00000087);   // this corresponds to Q(x) = x^7 + x^2 + x^1 + 1
    __m128i k = _mm_xor_si128(_mm_clmulepi64_si128(b, q, 0x01), _mm_slli_si128(b, 8)); // K = CLMUL(B[1], Q) + B[0]*x^64

    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1k0 = _mm_clmulepi64_si128(a, k, 0x01);
    __m128i a1k1 = _mm_clmulepi64_si128(a, k, 0x11);

    __m128i lower  = _mm_xor_si128(a0b0, a1k0);     // computes lower  = a0b0 + a1k0
    __m128i higher = _mm_xor_si128(a0b1, a1k1);     // computes higher = a0b1 + a1k1

    __m128i c01 = _mm_xor_si128(lower, _mm_slli_si128(higher, 8));  // c01 = lower + (higher << 64)
    __m128i c23 = _mm_srli_si128(higher, 8);                        // c02 = higher >> 64           (highest 64 bits are 0)

    // Step 2: Reduce
    __m128i x = _mm_clmulepi64_si128(c23, q, 0x00);
    c01 = _mm_xor_si128(x, c01);

    return c01;
}
```

## GHASH capable GFMUL implementation
Now that we have two functional versions of computing the multiplication of two GF($2^{128}$) elements, let's see how we can use it inside the GHASH calculation required for AES-GCM.

The NIST specification declares that their polynomials are stored exactly the opposite way round of how we store them. 
- We store $h(x) = x^3 + x^2 + x^1$ as follows: `__m128i h = _mm_set_epi32(0, 0, 0, 0x0000000e);` (0..000001110)
- NIST specifies it the other way round: `__m128i h = _mm_set_epi32(0x70000000, 0, 0, 0);` (01110000..0)
and that all GHASH operations happen on bit-reflected values.

Thus, to convert our currently used polynomials to abide to the NIST specification, we have to apply the bit reflection transformation: $R(..)$ [see the mathematical section](#bit-reflection).

Given a block of actual AES-GCM data $\big(a(x)\big)$ and the GHASH key $H$ $\big(b(x)\big)$, in order to compute their galois-field multiplication result $\big(d(x)\big)$, we would have to do:
```math
d(x) := R\Big(\text{GFMUL}\big(R(a(x)), R(b(x))\big)\Big)
```
where $\text{GFMUL}$ could be either our standard or the K-optimized implementation.

To avoid costly bit-reflection, we can use the CLMUL-identity which simply states:
```math
\text{CLMUL}\big(R_{128}(a(x)), R_{128}(b(x))\big) = R_{256}\big(\text{CLMUL}(a(x), b(x)) << 1\big) 
```
(under the assumption that $\text{CLMUL(a(x), b(x)) << 1}$ results in a polynomial of degree $255$, which is the maximum possible for that calculation and that $a(x)$ and $b(x)$ have coefficients $=0$ for all remaining terms if their degree is $\le 127$).

Important note: The following mathematical inscription of the operation might not be 100% correct and therefore contain some magic and inconsistencies.
It arouse from trying to understand the algorihm in "Table 9" in Intel Doc. B.
And this part is not actually complete but it might give hints on how to get to the final result.

Therefore, Applying the bit-reflection transformation on a polynomial, when simply writing $R(..)$, the chosen $n$ might dynamically expand to the needed resolution.
E.g.: Even though we write $R_{128}\Big(R_{256}\big(\text{CLMUL}(a(x), b(x)) << 1\big) \mod P(x)\Big)$, indicating that $R_{128}$ is applied on $P(x)$, applying $R_{129}$ fits our case better, especially since the degree of $P(x)$ is $128 \rightarrow n=129$.
Thus:
```math
\begin{align*}
d(x) &= R\Big(\text{GFMUL}\big(R(a(x)), R(b(x))\big)\Big) \\[6pt]
&= R_{128}\Big(\text{CLMUL}\big(R_{128}(a(x)), R_{128}(b(x))\big) \mod P(x)\Big) \\[6pt]
&\equiv R_{128}\Big(R_{256}\big(\text{CLMUL}(a(x), b(x)) << 1\big) \mod P(x)\Big) \\[6pt]
\end{align*}
```
Say:
```math
\begin{align*}
c^{*}(x) &= CLMUL(a(x), b(x)) << 1 \\[6pt]
&= c(x) << 1\\[6pt]
&= (c_{[3]}(x)x^{192} + c_{[2]}(x)x^{128} + c_{[1]}(x)x^{64} + c_{[0]}(x)) \cdot x \\[6pt]
&= c_{[3]}(x)x^{193} + c_{[2]}(x)x^{129} + c_{[1]}(x)x^{65} + c_{[0]}(x)x \\[6pt]
&= c^{*}_{[3]}(x)x^{192} + c^{*}_{[2]}(x)x^{128} + c^{*}_{[2]}(x)x^{64} + c^{*}_{[0]}(x) \\[6pt]
\end{align*}
```
where all $`c^{*}_{[i]}(x)`$ polynomials are of degree $\leq 63$, because the degree of $`c_{[3]}(x)`$ is $\leq$ 62$.\
Now, let's define: 
```math
A(x) := x^{255}c^{*}(x^{-1}) = x^{255}(\text{CLMUL}(a, b) << 1)(x^{-1})
```
and let $B$ be the remainder: 
```math
B(x) := A(x) \mod P(x)
```
Thus $A(x)$ can also be written in the following way for an unknown $L(x)$: 
```math
A(x) = P(x)\cdot L(x) + B(x)
```
It follows:
```math
B(x) = A(x) - P(x)\cdot L(x)
```
Let's insert this into the other equation from before:
```math
\begin{align*}
d(x) &= R_{128}\big(R_{255}(c^{*})(x^{-1} \mod P(x)\big) \\[6pt]
&\equiv R_{128}(B)(x) \\[6pt]
&\equiv x^{127}B(x^{-1}) \\[6pt]
&\equiv x^{127}(A - P\cdot L)(x^{-1}) \\[6pt]
&\equiv x^{127}A(x^{-1}) - x^{128}P(x^{-1})L(x^{-1}) \\[6pt]
&\equiv x^{127}x^{-255}c^{*}(x) - x^{128}P(x^{-1})L(x^{-1}) \\[6pt]
&\equiv x^{-128}c^{*}(x) - x^{128}(x^{-128} + x^{-7} + x^{-2} + x^{-1} + x^{-0})L(x^{-1}) \\[6pt]
&\equiv x^{-128}(c^{*}_{[3]}(x)x^{192} + c^{*}_{[2]}(x)x^{128} + c^{*}_{[1]}(x)x^{64} + c^{*}_{[0]}(x)) - (x^{0} + x^{121} + x^{126} + x^{127} + x^{128})L(x^{-1}) \\[6pt]
&\equiv (c^{*}_{[3]}(x)x^{64} + c^{*}_{[2]}(x)x^{0} + c^{*}_{[1]}(x)x^{-64} + c^{*}_{[0]}(x)x^{-128} - P_{r}(x)L(x{^-1}) \\[6pt]
\end{align*}
```
where the change from $x^{127}$ to $x^{128}$ for applying on $P(x)$ is some magic, and with $P_{r}(x) = x^{128} + x^{127} + x^{126} + x^{121} + 1$.
This is all the mathematics I was able to deduce.

One important note as for the memory representation:<br>
Say we store:
```math
\begin{align*}
c(x) &= c_{[3]}(x)x^{192} + c_{[2]}(x)x^{128} + c_{127}x^{127} + c_{126}x^{126} + ... + c_{64}x^{64} + c_{63}x^{63} + ... + c_{0}x^{0} \\[8pt]
c_{[1]}(x) &= c_{127}x^{127} + c_{126}x^{127} + ... c_{64}x^{64} \\[8pt]
c_{[0]}(x) &= c_{63}x^{63} + ... + c_{0}x^{0} \\[8pt]
\end{align*}
```
in memory as follows: `__m128i c = _mm_set_epi32(c3, c2, c1, c0);`.
Then after applying the bitshift to the left, **of** the original coefficients of `c0` corresponding to $c_{\[0\]}(x)$, the highest coefficient $c_{63}$ now resides in the second word, `c1`. $(c_{63}x^{63} \cdot x = c_{63}x^{64})$
Therefore, after the bitshift, `c0` does not correspond to $c_{\[0\]}(x)$ anymore.<br>
But, when looking at $`c^{*}(x) = c^{*}_{[3]}(x)x^{192} + c^{*}_{[2]}(x)x^{128} + c^{*}_{[1]}(x)x^{64} + c^{*}_{[0]}(x)`$, then $`c^{*}_{[0]}(x)`$ does correspond to `c0` again correctly.

-----

So much for the mathematics part.<br>
Our used implementation retrieved from Intel Doc.B with the small addition of shifting `C[3:0]` one to the left before applying the reduction (NECESSARY), can be seen in the following picture.
It does though not align perfectly with the maths, or at least I haven't seen how to get from one to the other yet.

The difference:
- We multiply $`Q'(x) := x^{127} + x^{126} + x^{121}`$ with `C'[0]` (being our $`c^{*}_{[0]}(x)`$), add it onto `C'[3:1]` and also add `C'[0]*x^64`. (This part might come from the $+1$ in $P_{r}(x)$).
- We repeat this step for the remaining part to then get our final result.

This complete sketch can be viewed here:
<img width="900" height="1246" alt="grafik" src="https://github.com/user-attachments/assets/f8684b06-6229-4757-a222-b9e2ebbffac7" />

Please again note the different notation.
Even though the picture uses $A'$ and $B'$, they are in fact our plain $a(x)$ and $b(x)$, and the result is our usable GHASH result.

### Code Implementation
See a full code implementation below using Intel AVX instruction set:
```c
_m128i gfmul_reversed_bl_opt(__m128i a, __m128i b){
    __m128i Q_r = _mm_set_epi32(0, 0, 0xc2000000, 0);   // Q' = Q_r = x^127+x^126+x^121

    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);

    __m128i mid = _mm_xor_si128(a0b1, a1b0);      // computes mid = A0B1 + A1B0

    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8));    // computes C[1:0] = A0B0 + (mid << x^64)
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8));    // computes C[3:2] = A1B1 + (mid >> x^64)

    // Step 1.1: Bitshift << 1
    __m128i tmp7,tmp8,tmp9;
    tmp7 = _mm_srli_epi32(c01, 31);
    tmp8 = _mm_srli_epi32(c23, 31);
    c01 = _mm_slli_epi32(c01, 1);
    c23 = _mm_slli_epi32(c23, 1);
    tmp9 = _mm_srli_si128(tmp7, 12);
    tmp8 = _mm_slli_si128(tmp8, 4);
    tmp7 = _mm_slli_si128(tmp7, 4);
    c01 = _mm_or_si128(c01, tmp7);
    c23 = _mm_or_si128(c23, tmp8);
    c23 = _mm_or_si128(c23, tmp9);


    // Step 2.1: Reduce C[0]
    __m128i x = _mm_clmulepi64_si128(c01, Q_r, 0x00);
    c23 = _mm_xor_si128(c23, _mm_srli_si128(x, 8));       // add higher half of x (X[1]) to lower part of C[3:2]
    c23 = _mm_xor_si128(c23, _mm_unpacklo_epi64(c01, ZERO));    // add only C[0] to lower part of C[3:2]  (zeroing out C[1]) (this corresponds to the x^0 part of Q_r)

    c01 = _mm_xor_si128(c01, _mm_slli_si128(x, 8));       // add lower half of x (X[0]) to upper part of C[1:0]

    // Step 2.2: Reduce C[1]
    x = _mm_clmulepi64_si128(c01, Q_r, 0x01);               // computes C[1] * Q = higher(C[1:0]) * lower(Q)
    c23 = _mm_xor_si128(c23, x);                          // add full X on C[3:2]
    c23 = _mm_xor_si128(c23, _mm_unpackhi_epi64(ZERO, c01)); // add C[1] on higher C[3:2] (zeroing out C[0]) (this corresponds to the x^0 part of Q_r)

    return c23;
}
```
