# Sparse Matrices

While Julia provides a data structure for sparse matrix computations, the
employed design decisions make it difficult to use this implementation for
computations over finite fields. This is mainly due to the fact that in
the Julia implementation, it is assumed that one can determine the zero and
one elements from the data type alone. However, a finite field data type
and especially the ones implemented in the package Nemo, generally also
depends on additional parameters, such as the characteristic of the field.

Since the algorithms underlying `ConleyDynamics.jl` only require basic row
and column operations, we decided to include a specialized sparse matrix
implementation.

## Sparse Matrix Format



