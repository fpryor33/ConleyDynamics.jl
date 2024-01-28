# ConleyDynamics.jl

*Conley index and multivector fields for Julia.*

A Julia package for studying combinatorial multivector fields using Conley theory

!!! note

    This package and the documentation is still very much work in progress!

## Package Features

- Implements data structures for Lefschetz complexes, in particular simplicial
  and cubical complexes.
- Classical Forman combinatorial vector fields and multivector fields are
  supported.
- Computation of connection matrices and Conley-Morse graphs.
- Basic homology algorithms over finite fields, including persistent homology
  and relative homology.
- Includes a sparse matrix implementation geared towards computations over
  finite fields.

We can also include Latex, as in ``u_t = -\Delta(\epsilon^2 \Delta u + f(u))``.

```@docs
LefschetzComplex
```

Here is another entry:

```@docs
ConleyMorseCM
```

And another one:

```@docs
SparseMatrix
```



