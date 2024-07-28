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

## Definition of a Lefschetz Complex

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

At first glance the above definition can seem daunting. However,
it is based on a straightforward geometric idea. A Lefschetz 
complex is a structure that is built from elementary building
blocks called *cells*. Each cell has a dimension associated with
it, and it is topologically an open ball of this dimension. Thus,
cells of dimension zero are points, also called *vertices*. Cells
of dimension one are open curve segments, which we call *edges*,
and two-dimensional cells are called *faces* and take the form
of open two-dimensional membranes.

The incidence coefficient map encodes how these cells are glued
together to form the Lefschetz complex ``X``. In order to shed
more light on this gluing process, counsider the *boundary map*
``\partial`` which is defined on cells via

```math
   \partial x = \sum_{y \in X} \kappa(x,y) y \; .
```

This map sends a cell ``x`` of dimension ``k`` to a specific 
linear combination of cells of dimension ``k-1``, called the
*boundary* of ``x``. By using ideas from linear algebra, the
boundary map can be extended to map a general linear combination
of ``k``-dimensional cells to the corresponding linear combination
of the separate boundaries. For example, if one chooses the field
``F = \mathbb{Q}`` of rationals, one has ``\partial  (x_1 - 2x_2)
= \partial x_1 - 2 \partial x_2``. Notice that using this extended
definition of the boundary map, one can rewrite the summation
condition in the definition of a Lefschetz complex in the
equivalent form

```math
   \partial( \partial x) = 0
   \quad\text{ for all cells }\quad
   x \in X \; .
```

In other words, the boundary of any cell is itself boundaryless.

With the help of the boundary map, one can frequently infer the
overall geometric structure of a Lefschetz complex ``X``. For this,
think of a Lefschetz complex as being build *from the ground up* in
the following way. First, start by putting down all vertices of ``X``
at different locations in some ambient space. Since the boundary of
each one-dimensional cell is made up of a linear combination of
vertices, one can then add a curve segment for each one-dimensional
cell, which connects the vertices in its boundary. Note that in the
general version of a Lefschetz complex it is possible that an edge
has only one vertex in its boundary, or maybe even none, and in these
cases the edge is either only connected to the one boundary vertex,
or it is an open curve segment connected to no vertex at all,
respectively. Continue in this fashion to add two-dimensional
faces to fill in the space between the edges in its boundary,
and so on for higher dimensions. Needless to say, in the case
of a general complicated Lefschetz complex this procedure is
of limited use, since the boundary of a cell can be an arbitrary
linear combination of cells, with coefficients that can be any
nonzero numbers in the field ``F``. Yet, in many simple cases
the above intuition is sufficient.

## Lefschetz complex data structure

In `ConleyDynamics.jl` a Lefschetz complex is implemented as the
following specific composite data type:

```@docs; canonical=false
LefschetzComplex
```

The fields of this struct relate to the mathematical definition
of a Lefschetz complex ``X`` in the following way:

- The integer `ncells` gives the total number of cells in ``X``.
  Internally, these cells are numbered by integers ranging from 1
  to `ncells`.
- The vector `dimensions` is a `Vector{Int}` and collects the 
  dimensions of the cells. In other words, the cell which is indexed
  by `k` has dimension `dimensions[k]`.
- The integer `dim` describes the overall dimension of the Lefschetz
  complex, which is the largest dimension of a cell.
- The incidence coefficient map ``\kappa`` is encoded in the sparse
  matrix `boundary`. This matrix is a square matrix with `ncells` 
  rows and columns. The ``k``-th column contains the incidence
  coefficients ``\kappa(k,\cdot)`` in the sense that ``\kappa(k,m)``
  equals the entry in row ``m`` and column ``k``. Since for most
  Lefschetz complexes the majority of the incidence coefficients
  is zero, the matrix is represented using the sparse format
  [`SparseMatrix`](@ref), which is described in more detail
  in [Sparse Matrices](@ref).
- While the internal representation of cells as integers is 
  computationally convenient, it does make interpreting the
  results more cumbersome. Each Lefschetz complex therefore has
  to have string labels assigned to each cell as well. These are
  contained in `labels::Vector{String}`, where `labels[k]` gives
  the label of cell `k`.
- In order to easily determine the integer index for a cell with
  a specific label, the field `indices` contains a dictionary
  of type `Dict{String,Int}` which maps labels to indices. For 
  example, if a cell has the label `"124.010"`, then the associated
  integer index is given by `indices["124.010"]`.

A object of type `LefschetzComplex` is created by passing the
field items in the order given in [`LefschetzComplex`](@ref).
Consider for example the Lefschetz complex from Figure 4
in [mrozek:wanner:p21a](@cite), see also the left complex in the
next image. This complex consists of six cells, and we initialize
the vector of labels, the cell index dictionary, and the cell 
dimensions via the commands

```julia
ncL = 6
labelsL  = Vector{String}(["A","B","a","b","c","alpha"])
indicesL = Dict{String,Int}([(labelsL[k],k) for k in 1:length(labelsL)])
cdimsL   = [0, 0, 1, 1, 1, 2]
```

The boundary matrix can then be defined using

```julia
bndmatrixL = zeros(Int, ncL, ncL)
bndmatrixL[[1,2],3] = [1; 1]     # a
bndmatrixL[[1,2],4] = [1; 1]     # b
bndmatrixL[[1,2],5] = [1; 1]     # c
bndmatrixL[[3,4],6] = [1; 1]     # alpha
bndsparseL = sparse_from_full(bndmatrixL, p=2)
```

Notice that we first create the matrix as a regular integer 
matrix, and then use the function [`sparse_from_full`](@ref) 
to turn it into sparse format over the field ``GF(2)`` with
characteristic `p = 2`. Finally, the Lefschetz complex is
created using

```julia
lcL = LefschetzComplex(ncL, 2, bndsparseL, labelsL, indicesL, cdimsL)
```

![Two sample Lefschetz complexes](img/lefschetzex1.png)

Lefschetz complexes do not always have to contain cells of
all dimensions. For example, the Lefschetz complex shown on the
right side of the figure has no vertices, and it can be created
using the commands

```julia
ncR = 4
labelsR  = Vector{String}(["a","b","c","alpha"])
indicesR = Dict{String,Int}([(labelsR[k],k) for k in 1:length(labelsR)])
cdimsR   = [1, 1, 1, 2]
bndmatrixR = zeros(Int, ncR, ncR)
bndmatrixR[[1,2,3],4] = [1; 1; 1]     # alpha
bndsparseR = sparse_from_full(bndmatrixR, p=2)
lcR = LefschetzComplex(ncR, 2, bndsparseR, labelsR, indicesR, cdimsR)
```

## Basic concepts in Lefschetz complexes



## Simplicial complexes


## Cubical complexes





[dlotko:etal:11a](@cite)
[lipinski:etal:23a](@cite)





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

