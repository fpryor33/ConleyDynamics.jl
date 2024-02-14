# Examples

```@meta
DocTestSetup = quote
    push!(LOAD_PATH,"../../../src/")
    using ConleyDynamics
end
```

In the following we discuss a number of connection matrix examples.

## Critical flow on a simplex

```@docs
example_critical_simplex(::Int)
```

## Combinatorial flows on a cylinder and a Moebius strip

```@docs
example_moebius()
```

## Nonunique connection matrices

![An example with nonunique connection matrices](img/multiconn.png)

```@docs
example_nonunique()
```

## Examples from MW-2023

The following examples are taken from [mrozek:wanner:p21a](@cite).

![Four sample combinatorial vector fields](img/connectionex.png)

```@docs
example_MW_fig02()
```

```@docs
example_MW_fig01()
example_MW_fig03()
example_MW_fig04()
example_MW_fig11()
```

## Examples from BKMW-2020

The following examples are taken from [batko:etal:20a](@cite).

```@docs
example_BKMW20_fig1()
example_BKMW20_fig3()
```

## Examples References

See the [full bibliography](@ref References) for a complete list
of references cited throughout this documentation. This section cites
the following references:

```@bibliography
Pages = ["examples.md"]
Canonical = false
```

