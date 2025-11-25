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


# GFMUL Implementation
\label{label}
----
This document is structured into the following key points:
0. [Mathematical notations](#mathematical-notations)
1. [Standard GMFUL implementation](#standard-gfmul-implementation)
2. [K-Optimization](#k-optimization)
3. [GHASH capable GFMUL implementation](#ghash-capable-gfmul-implementation)

----

## Mathematical Notations
GFMUL stands for Galois-Field multiplication. 
For AES-GCM this refers to the GF($2^128$) field.
Elements of this field can be expressed as polynoms with coefficients in GF($2^1$)
```math
a \in \text{GF}(2^128) = \sum_{i=0}^{127}a_i \cdot x^i, \text{with} a_i \in \text{GF}(2^1)
```

## Standard GFMUL implementation

## K-Optimization


## GHASH capable GFMUL implementation
