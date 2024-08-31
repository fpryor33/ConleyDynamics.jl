# Examples

```@meta
DocTestSetup = quote
    push!(LOAD_PATH,"../../../src/")
    using ConleyDynamics
end
```

In order to illustrate the basic functionality of `ConleyDynamics.jl`, 
this section collects a number of examples. Many of these are taken 
from the papers [batko:etal:20a](@cite) and [mrozek:wanner:p21a](@cite),
and they consider both Forman vector fields and general multivector
fields on a variety of Lefschetz complexes. Each example has its own
associated function, so that users can quickly create examples
on their own by taking the respective source files as templates.

## Critical Flow on a Simplex

The first example considers the arguably simplest situation of a
Forman vector field on a simplicial complex. The simplicial complex
``X`` is given by a single simplex of dimension ``n``, together
with all its faces, while the Forman vector field on ``X`` contains
only singletons. In other words, every simplex in the complex
is a critical cell. Thus, this combinatorial dynamical system has
one equilibrium of index ``n``, and ``n+1`` stable equilibria. In
addition, there are ``2^{n+1} - n - 3`` additional stationary states
whose indices lie strictly between ``0`` and ``n``, as well as a 
wealth of algebraically induced heteroclinic orbits. All of these
can be found by using the connection matrix for the problem, as
outlined in the following description for the function
[`example_critical_simplex`](@ref).

```@docs; canonical=false
example_critical_simplex(::Int)
```

## Flow on a Moebius Strip

```@docs; canonical=false
example_moebius
```

## Nonunique Connection Matrices

![An example with nonunique connection matrices](img/multiconn.png)

```@docs; canonical=false
example_nonunique()
```

## Further Connection Matrix Examples

The following examples are taken from [mrozek:wanner:p21a](@cite).

![Four sample combinatorial vector fields](img/connectionex.png)

```@docs; canonical=false
example_MW_fig02()
```

```@docs; canonical=false
example_MW_fig01()
```

```@docs; canonical=false
example_MW_fig03()
```

```@docs; canonical=false
example_MW_fig04()
```

```@docs; canonical=false
example_MW_fig11()
```

## Forman Vector Field Examples

The following examples are taken from [batko:etal:20a](@cite).

```@docs; canonical=false
example_BKMW20_fig1()
```


```@docs; canonical=false
example_BKMW20_fig3()
```



## [References](@id refexamples)

See the [full bibliography](@ref References) for a complete list
of references cited throughout this documentation. This section cites
the following references:

```@bibliography
Pages = ["examples.md"]
Canonical = false
```

