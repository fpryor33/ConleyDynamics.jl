# ConleyDynamics.jl

*Conley index and multivector fields for Julia.*

A Julia package for studying combinatorial multivector fields using Conley theory

!!! note

    This package and the documentation is still very much work in progress!

## Package Features

- Data structures for Lefschetz complexes, in particular simplicial and
  cubical complexes.
- Classical Forman combinatorial vector fields and multivector fields are
  supported.
- Computation of connection matrices and Conley-Morse graphs.
- Basic homology algorithms over finite fields, including persistent homology
  and relative homology.
- Algorithms rely on a built-in sparse matrix implementation geared towards
  computations over finite fields.

The [Package Guide](@ref) provides a tutorial explaining how to get started
using ConleyDynamics, and further examples can be found on the [Examples](@ref)
page.

## Manual Outline

```@contents
Pages = [
    "man/guide.md",
    "man/examples.md",
    "man/datastruct.md",
    "man/connections.md",
    "man/homology.md",
    "man/sparse.md",
]
Depth = 1
```

