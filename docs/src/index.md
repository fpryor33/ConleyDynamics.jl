# ConleyDynamics.jl

*Conley index and multivector fields for Julia.*

## Introduction

[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl) is a
Julia package for studying combinatorial multivector fields using Conley
theory. The multivector fields can be studied on arbitrary Lefschetz
complexes, which include both simplicial and cubical complexes as
important special cases. The concept of combinatorial multivector field
generalizes Forman vector fields, which were originally introduced to study
Morse theory in a discrete combinatorial setting.

!!! note
    This documentation is also available in PDF format: [ConleyDynamics.pdf](ConleyDynamics.pdf).

## Features

- Data structures for Lefschetz complexes, in particular simplicial and
  cubical complexes.
- Classical Forman combinatorial vector fields and multivector fields are
  supported.
- Computation of Conley indices, connection matrices, and Conley-Morse graphs.
- Basic homology algorithms over finite fields and the rationals, including
  persistent homology and relative homology.
- Algorithms rely on a built-in sparse matrix implementation which is geared
  towards computations over finite fields and the rationals.

## Installation

To use [ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl)
please install Julia 1.10 or higher. See <https://julialang.org/downloads/>
for instructions on how to obtain Julia for your system.

At the Julia prompt simply type

```
julia> using Pkg; Pkg.add("ConleyDynamics")
```

After Julia has finished downloading and precompiling the package and
all of its dependencies, you can start using it by typing

```
julia> using ConleyDynamics
```

## Manual Outline

The [Tutorial](@ref) briefly explains how to get started with
[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl).
More details, including on the underlying mathematics, are provided in
the following three sections, which cover Lefschetz complexes, homology,
and Conley theory including connection matrices. After a discussion of all
included examples in the [Examples](@ref) section, the manual concludes
with a description of the sparse matrix format underlying the package.

```@contents
Pages = [
    "man/tutorial.md",
    "man/lefschetz.md",
    "man/homology.md",
    "man/conley.md",
    "man/examples.md",
    "man/sparse.md",
]
Depth = 2
```

