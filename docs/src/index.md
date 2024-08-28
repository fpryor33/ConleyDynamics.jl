# ConleyDynamics.jl

*Conley index and multivector fields for Julia.*

## Introduction

`ConleyDynamics.jl` is a Julia package for studying combinatorial multivector
fields using Conley theory. The multivector fields can be studied on arbitrary
Lefschetz complexes, which include both simplicial and cubical complexes
as important special cases. The concept of combinatorial multivector field
generalizes Forman vector fields, which were originally introduced to study
Morse theory in a discrete combinatorial setting.

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

To use `ConleyDynamics.jl` please install Julia 1.10 or higher. See
<https://julialang.org/downloads/> for instructions on
how to obtain Julia for your system.

At the Julia prompt simply type

```
julia> using Pkg; Pkg.add("ConleyDynamics")
```

## Manual Outline

The [Tutorial](@ref) briefly explains how to get started with
`ConleyDynamics.jl`. More details, including on the underlying mathematics,
are provided in the following four sections, which cover Lefschetz complexes,
homology, Conley theory including connection matrices, and sparse matrices,
respectively. Finally, a discussion of all included examples can be found in the
[Examples](@ref) section.

```@contents
Pages = [
    "man/tutorial.md",
    "man/lefschetz.md",
    "man/homology.md",
    "man/conley.md",
    "man/sparse.md",
    "man/examples.md",
]
Depth = 2
```

