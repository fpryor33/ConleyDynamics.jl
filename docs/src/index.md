# ConleyDynamics.jl

*Conley index and multivector fields for Julia.*

## Introduction

ConleyDynamics.jl is a Julia package for studying combinatorial multivector
fields using Conley theory. The multivector fields can be studied on arbitrary
Lefschetz complexes, which include both simplicial and cubical complexes
as important special cases. The concept of combinatorial multivector field
generalizes Forman vector fields, which were originally introduced to study
Morse theory in a discrete combinatorial setting.

!!! note

    This package and the documentation is still very much work in progress!

## Features

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

## Installation

To use ConleyDynamics please install Julia 1.10 or higher. See
<https://julialang.org/downloads/> for instructions on
how to obtain Julia for your system.

At the Julia prompt simply type

```
julia> using Pkg; Pkg.add("ConleyDynamics")
```

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

