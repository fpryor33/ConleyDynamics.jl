# Conley Theory

Here we need a more detailed description of Conley theory, and in particular
connection matrices.  In particular, this should discuss the various field types
that can be used, as well as all the entries in the [`ConleyMorseCM`](@ref) data
structure.

[kaczynski:etal:16a](@cite)
[mrozek:wanner:21a](@cite)


## Intervals in a planar Morse decomposition

We briefly return to one of the examples in the tutorial. More precisely, we
consider the planar ordinary differential equation given by

```math
   \begin{array}{rcl}
     \dot{x}_1 & = & x_1 \left( 1 - x_1^2 - 3 x_2^2 \right) \\[1ex]
     \dot{x}_2 & = & x_2 \left( 1 - 3 x_1^2 - x_2^2 \right)
   \end{array}
```

The dynamics of this system is characterized by the existence of a global
attractor in the shape of a closed disk. Inside the attractor, there are
nine different Morse sets:

- The origin is an equilibrium of index 2, i.e., it is an unstable
  stationary state with a two-dimensional unstable manifold.
- The four points ``(\pm 1/2, \pm 1/2)`` are unstable equilibria of
  index 1, i.e., with a one-dimensional unstable manifold.
- Finally, the four points ``(\pm 1, 0)`` and ``(0, \pm 1)`` are
  asymptotically stable stationary states.

We saw in the tutorial that the Morse decomposition of this system
can easily be found using `ConleyDynamics.jl`, as well as the
associated connection matrix. Yet, in certain situations one might
only be interested in part of the dynamics on the attractor. Moreover,
while the Morse sets describe the recurrent part of the dynamics, they
do not provide information on the geometry of the connecting sets between
the Morse sets. In the following, we illustrate how this can be analyzed
further.

The right-hand side of the above vector field can be implemented using
the Julia function

```@example Cmorseinterval
using ..ConleyDynamics # hide
function planarvf(x::Vector{Float64})
    #
    # Sample planar vector field with nontrivial Morse decomposition
    #
    x1, x2 = x
    y1 = x1 * (1.0 - x1*x1 - 3.0*x2*x2)
    y2 = x2 * (1.0 - 3.0*x1*x1 - x2*x2)
    return [y1, y2]
end
```

To analyze the resulting global dynamical behavior, we first create a
simplicial mesh covering the square ``[-6/5, 6/5]^2`` using the commands

```@example Cmorseinterval
lc, coords = create_simplicial_delaunay(300, 300, 5, 50);
coordsN = convert_planar_coordinates(coords,[-1.2,-1.2], [1.2,1.2]);
lc.ncells
```

The integer in the output gives the number of cells in the created Lefschetz
complex ``X``. Note that we are using a Delaunay triangulation over an initial
box of size ``300 \times 300``, where the target triangle size is about 5.
This box is then rescaled to cover the above square. We can then create a 
multivector field on the simplicial complex `lc` and find its Morse
decomposition using the commands

```@example Cmorseinterval
mvf = create_planar_mvf(lc, coordsN, planarvf);
morsedecomp = morse_sets(lc, mvf);
length(morsedecomp)
```

As expected, `ConleyDynamics.jl` finds exactly nine Morse sets. Their 
Conley indices can be computed and stored in a `Vector{Vector{Int}}`
using the command

```@example Cmorseinterval
conleyindices = [conley_index(lc, mset) for mset in morsedecomp]
```

These Conley indices correspond to the dynamical behavior near the
equilibrium solutions described above.

Suppose now that rather than finding the connection matrix for the
complete Morse decomposition, we would only like to consider a 
part of it. This can be done as long as we restrict our attention
to an *interval in the Morse decomposition*. Such an interval
``\mathcal{I}`` can be created from a selection ``\mathcal{S}``
of the Morse sets in the following way:

- In addition to the Morse sets in ``\mathcal{S}``, the interval
  ``\mathcal{I}`` contains all Morse sets that lie between two
  Morse sets in ``\mathcal{S}`` with respect to the poset order
  underlying the Morse decomposition. Recall that this poset order
  can be computed via [`morse_sets`](@ref) by activating the extra
  return object `hasse`, which describes the Hasse diagram of the
  poset.

With every interval ``\mathcal{I}`` of the Morse decomposition one
can assign a smallest isolated invariant set ``X_{\mathcal{I}}
\subset X`` which describes the complete dynamics within and
between the Morse sets in ``\mathcal{I}``. In fact, in can be
characterized as follows:

- The set ``X_{\mathcal{I}}`` consists of all cells in the 
  underlying Lefschetz complex ``X`` through which one can find
  a solution which originates in one Morse set of ``\mathcal{I}``
  and ends in another Morse set of ``\mathcal{I}``, where the
  two involved Morse sets can be the same. In other words, one
  needs to combine the interval Morse sets with all *connecting
  orbits* between them.

The two above steps can be performed in `ConleyDynamics.jl` using
the function [`morse_interval`](@ref).

In our example, we consider two intervals. The first interval
consists of the five Morse sets corresponding to all unstable
equilibrium solutions, while the second one considers the
four index 1 and the four stable stationary states. The
associated isolated invariant sets for these two intervals
can be computed as follows:

```@example Cmorseinterval
subset1 = findall(x -> x[2]+x[3]>0, conleyindices);
subset2 = findall(x -> x[1]+x[2]>0, conleyindices);
lcsub1 = morse_interval(lc, mvf, morsedecomp[subset1]);
lcsub2 = morse_interval(lc, mvf, morsedecomp[subset2]);
[length(subset1), length(subset2), length(lcsub1), length(lcsub2)]
```

The output shows that we have in fact extracted five and eight
Morse sets, respectively. It also shows that the Lefschetz
complexes corresponding to these two isolated invariant sets
are much smaller than ``X``.

So far, we have just determined the collections of cells that 
correspond to the two isolated invariant sets for these intervals.
We can now restrict the combinatorial dynamics to these subsets. 
Note that since they are both isolated invariant sets, they are
locally closed in ``X``, and therefore the restrictions provide
us with two new Lefschetz complexes `lcr1` and `lcr2`, along with
induced multivector fields `mvfr1` and `mvfr2`, respectively. In
`ConleyDynamics.jl`, this is achieved using the commands

```@example Cmorseinterval
lcr1, mvfr1 = restrict_dynamics(lc, mvf, lcsub1);
lcr2, mvfr2 = restrict_dynamics(lc, mvf, lcsub2);
[lcr1.ncells, lcr2.ncells]
```

It is now easy to find the connection matrices for these two
intervals. The first connection matrix is given by

```@example Cmorseinterval
cmr1 = connection_matrix(lcr1, mvfr1);
cmr1.conley
```

```@example Cmorseinterval
full_from_sparse(cmr1.matrix)
```

It clearly shows that the unstable index 2 Morse set
has connecting orbits to every one of the four index 1 
equilibria. Similarly, the second connection matrix can
be determined as

```@example Cmorseinterval
cmr2 = connection_matrix(lcr2, mvfr2);
cmr2.conley
```

```@example Cmorseinterval
full_from_sparse(cmr2.matrix)
```

In this case, every index 1 equilibrium is connected 
two its two neighboring stable stationary states via
heteroclinics that are detected by the connection matrix.

The Lefschetz complexes associated with the two Morse
decomposition intervals can also be visualized in
`ConleyDynamics.jl`. For this, recall that the function
[`plot_planar_simplicial_morse`](@ref) can plot an underlying
simplicial complex together with any collection of cell 
subsets. For our purposes, we use the following commands:

```julia
show1 = [[lcr1.labels]; cmr1.morse];
show2 = [[lcr2.labels]; cmr2.morse];
fname1 = "/Users/wanner/Desktop/invariantinterval2d1.png"
fname2 = "/Users/wanner/Desktop/invariantinterval2d2.png"
plot_planar_simplicial_morse(lc, coordsN, fname1, show1, vfac=1.1, hfac=2.0)
plot_planar_simplicial_morse(lc, coordsN, fname2, show2, vfac=1.1, hfac=2.0)
```

The variable `show1` collects not only the Morse sets that
are part of the first connection matrix `cmr1`, but also the
support of the Lefschetz complex `lcr1`. This support is accessed
via `[lcr1.labels]`, and we add it as a first vector of cells
in `show1`. Similarly, we determine the support of the second 
isolated invariant set, together with the Morse sets of `cmr2`.
The remaining four commands create two images.

![Interval support for the first interval](img/invariantinterval2d1.png)

The first image shows the five Morse sets surrounding the stationary
states at the origin and at ``(\pm 1/2, \pm 1/2)``. In addition, it
highlights the support of the isolated invariant set associated with
this Morse decomposition interval. One can clearly see rough outer
approxomations for the four heteroclinics which start at the origin
and end at the index 1 equilibria. These approximations are necessarily
coarse, since we are not working with a very fine triangulation.

![Interval support for the second interval](img/invariantinterval2d2.png)

Finally, the second image depicts the eight Morse sets enclosing
the index 1 and the stable stationary states. It also shows the
support of the Lefschetz complex `lcr2` which is associated with this
Morse decomposition interval. In this case, it covers eight different
heteroclinic orbits, which are in fact better approximated than the
four in the previous image.

## Analyzing a Planar System using Cubical Complexes

Our next example illustrates how `ConleyDynamics.jl` can be used to analyze
the global dynamics of a planar ordinary differential equations. For this,
consider the planar system

```math
   \begin{array}{rcl}
     \dot{x}_1 & = &  x_2 - x_1 \left( x_1^2 + x_2^2 - 4 \right)
       \left( x_1^2 + x_2^2 - 1 \right) \\[1ex]
     \dot{x}_2 & = & -x_1 - x_2 \left( x_1^2 + x_2^2 - 4 \right)
       \left( x_1^2 + x_2^2 - 1 \right)
   \end{array}
```

This system has already been considered in [mrozek:etal:22a](@cite).
The right-hand side of this vector field can be implemented using
the Julia function

```@example Ccircle
using ..ConleyDynamics # hide
function circlevf(x::Vector{Float64})
    #
    # Sample vector field with nontrivial Morse decomposition
    #
    x1, x2 = x
    c0 = x1*x1 + x2*x2
    c1 = (c0 - 4.0) * (c0 - 1.0)
    y1 =  x2 - x1 * c1
    y2 = -x1 - x2 * c1
    return [y1, y2]
end
```

To analyze the global dynamics of this vector field, we first create
a cubical complex covering the square ``[-3, 3]^2`` using the commands

```@example Ccircle
n = 51
lc, coords = create_cubical_rectangle(n,n,p=2);
coordsN = convert_planar_coordinates(coords,[-3.0,-3.0],[3.0,3.0]);
lc.ncells
```

As the last result shows, this gives a Lefschetz complex with 10609
cells. The multivector field can be generated using

```@example Ccircle
mvf = create_planar_mvf(lc, coordsN, circlevf);
length(mvf)
```

This multivector field consists of 2437 multivectors. Finally, the 
connection matrix can be determined using the command

```@example Ccircle
cm = connection_matrix(lc, mvf);
cm.conley
```

Therefore, the above planar system has three isolated invariant sets.
One has the Conley index of a stable equilibrium, while the other two
have that of a stable and an unstable periodic orbit. The columns of
the connection matrix correspond to these invariant sets as follows

```@example Ccircle
cm.poset
```

The connection matrix itself is given by

```@example Ccircle
full_from_sparse(cm.matrix)
```

This implies that there are connecting orbits from the unstable 
periodic orbit to both the stable equilibrium, and the stable
periodic orbit. To visualize these Morse sets, we employ the 
commands

```julia
fname = "cubicalcircles.pdf"
plot_planar_cubical_morse(lc, fname, cm.morse, pv=true)
```

![Morse sets of the planar circles vector field](img/cubicalcircles.png)

In the above example we used the original fixed cubical grid,
which is just a scaled version of the grid on the integer
lattice. It is also possible to work with a *randomized
grid*, in which the coordinates of the vertices are randomly
perturbed. This can be achieved with the following commands:

```julia
nR = 75
lcR, coordsR = create_cubical_rectangle(nR,nR,p=2,randomize=0.33);
coordsRN = convert_planar_coordinates(coordsR,[-3.0,-3.0],[3.0,3.0]);
mvfR = create_planar_mvf(lcR, coordsRN, circlevf);
cmR = connection_matrix(lcR, mvfR);
fnameR = "cubicalcirclesR.pdf"
plot_planar_cubical_morse(lcR, coordsRN, fnameR, cmR.morse, pv=true, vfac=1.1, hfac=2.0)
```

![Morse sets of the planar circles vector field via randomized cubes](img/cubicalcirclesR.png)

To constrast the above example with the use of a Delaunay triangulation,
we reanalyze the vector field in the following way:

```julia
lc2, coords2 = create_simplicial_delaunay(400, 400, 10, 30, p=2)
coords2N = convert_planar_coordinates(coords2,[-3.0,-3.0], [3.0,3.0])
mvf2 = create_planar_mvf(lc2, coords2N, circlevf)
cm2 = connection_matrix(lc2, mvf2)

fname2 = "cubicalcircles2.pdf"
plot_planar_simplicial_morse(lc2, coords2N, fname2, cm2.morse, pv=true)
```

In this case, the Morse sets can be visualized as in the figure.

![Morse sets of the planar circles vector field via Delaunay](img/cubicalcircles2.png)

Notice that we can also show the individual multivectors in more detail.
For the above example, we can plot all multivectors of the multivector
field `mvf2` which consist of at least 10 cells using the commands

```julia
mv_indices = findall(x -> (length(x)>9), mvf2)
large_mv = mvf2[mv_indices]

fname3 = "cubicalcircles3.pdf"
plot_planar_simplicial_morse(lc2, coords2N, fname3, large_mv, pv=true)
```

![Large multivectors in the Delaunay multivector field](img/cubicalcircles3.png)

Note that in this example, there are only 20 large multivectors.

## Conley Theory References

See the [full bibliography](@ref References) for a complete list
of references cited throughout this documentation. This section cites
the following references:

```@bibliography
Pages = ["conley.md"]
Canonical = false
```

