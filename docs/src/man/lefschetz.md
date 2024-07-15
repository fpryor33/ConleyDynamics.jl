# Lefschetz Complexes

The fundamental structure underlying the functionality of `ConleyDynamics.jl` is
a *Lefschetz complex*. It provides us with the basic model of phase space for
combinatorial dynamics. In view of the combinatorial, and therefore discrete,
character of the dynamical behavior, a Lefschetz complex is not a typical phase
space in the sense of classical dynamics. While the latter one is usually a
Euclidean space, a Lefschetz complex is basically a combinatorial model of it.
In the following fairly mathematical discussion, we provide its precise
mathematical definition, and explain how it can be created and modified within
the package. We also discuss two important special cases, namely *simplicial
complexes* and *cubical complexes*.

## Mathematical Definition of a Lefschetz Complex

The original definition of a Lefschetz complex can be found in
[lefschetz:42a](@cite), where it was simply referred to as a *complex*.

!!! tip "Definition: Lefschetz complex"
    Let ``F`` denote an arbitrary field. Then a pair ``(X,\kappa)``
    is called a *Lefschetz complex* over ``F`` if
    ``X = (X_k)_{k \in \mathbb{N}_0}`` is a finite set with
    ``\mathbb{N}_0``-gradation, and ``\kappa : X \times X \to F``
    is a mapping such that
    ```math
       \kappa(x,y) \neq 0
       \quad\mathrm{ implies }\quad
       x \in X_k
       \quad\mathrm{ and }\quad y \in X_{k-1},
    ```
    and such that for any ``x,y \in X`` one has
    ```math
       \sum_{z \in X} \kappa(x,z) \kappa(z,y) = 0 \; .
    ```
    The elements of ``X`` are referred to as *cells*, the
    value ``\kappa(x,y)`` is called the *incidence coefficient*
    of the cells ``x`` and ``y``, and the map ``\kappa`` is the
    *incidence coefficient map*. In addition, one defines the
    *dimension* of a cell ``x\in X_k`` as ``k``, and denotes it
    by ``k = \dim x``. Whenever the incidence coefficient map
    is clear from context, we often just refer to ``X`` as the
    *Lefschetz complex*.



[dlotko:etal:11a](@cite)
[lipinski:etal:23a](@cite)
[mrozek:wanner:p21a](@cite)




Here we need a more detailed description of Lefschetz complexes.
In particular, this should discuss the various field types that
can be used, as well as all the entries in the [`LefschetzComplex`](@ref)
data structure.



## Lefschetz Complexes References

See the [full bibliography](@ref References) for a complete list
of references cited throughout this documentation. This section cites
the following references:

```@bibliography
Pages = ["lefschetz.md"]
Canonical = false
```

