# Homology

Conley's theory for the qualitative study of dynamical systems is based
on fundamental concepts from algebraic topology. One of these is homology,
which studies the topological properties of spaces. As part of ConleyDynamics
a number of homology methods are included. These are described in detail 
below. More detailed information on the discussed topics can be found
in [munkres:84a](@citet*).

## Homology of Lefschetz Complexes

As we saw earlier, a Lefschetz complex is a collection of cells which
associated nonnegative dimensions, together with a boundary map.

```math
   \ldots \stackrel{\partial_{k+2}}{\longrightarrow}
   C_{k+1}(L) \stackrel{\partial_{k+1}}{\longrightarrow}
   C_{k}(L) \stackrel{\partial_{k}}{\longrightarrow}
   C_{k-1}(L) \stackrel{\partial_{k-1}}{\longrightarrow} \ldots
   \stackrel{\partial_{1}}{\longrightarrow}
   C_0(L) \longrightarrow 0
```

The boundary map ``\partial_{k}`` is a linear map from the chain group
``C_k(L)`` to the chain group ``C_{k-1}(L)``. Since ConleyDynamics uses
only field coefficients, these chain groups are in fact vector spaces. As
such, any linear map induces two important subspaces, which in the context
of algebraic topology are called as follows:

- The elements of the subspace ``Z_k(L) \subset C_k(L) = \mathrm{ker}\;
  \partial_k`` are called the _k-cycles_ of ``L``.
- The elements of the subspace ``B_k(L) \subset C_k(L) = \mathrm{im}\;
  \partial_{k+1}`` are called the _k-boundaries_ of ``L``.

Recall that the fundamental property of the boundary map is the identity
``\partial^2 = 0``, i.e., its square vanishes. Using the dimension-induced
gradation mentioned above, this can be reformulated as the equation
``\partial_{k} \circ \partial_{k+1} = 0``. This immmediately implies
the subspace inclusion ``B_k(L) \subset Z_k(L)``, and we can therefore
define the quotient space

```math
   H_k(L) =
   \mathrm{ker}\;\partial_k / \mathrm{im}\;\partial_{k+1}
```

This vector space is called the _k-th homology group_ of the Lefschetz 
complex ``L``.

## Relative Homology


## Persistent Homology of Filtrations


## Tutorial References

See the [full bibliography](@ref References) for a complete list
of references cited throughout this documentation. This section cites
the following references:

```@bibliography
Pages = ["homology.md"]
Canonical = false
```

