---
title: 'ConleyDynamics.jl: A Julia package for multivector dynamics on Lefschetz complexes'
tags:
  - Julia
  - dynamics
  - combinatorics
  - topology
  - mathematics
authors:
  - name: Thomas Wanner
    orcid: 0000-0003-3294-0366
    affiliation: 1
affiliations:
 - name: George Mason University, USA
   index: 1
date: 2 December 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Combinatorial topological dynamics is concerned with the qualitative study of
dynamical behavior on discrete combinatorial structures. It was originally developed
in the context of combinatorial vector fields [@forman:98a; @forman:98b] on simplicial
complexes, and has since been extended to the more general concept of combinatorial
multivector fields [@mrozek:17a; @lipinski:etal:23a] on arbitrary Lefschetz complexes
[@lefschetz:42a]. For such systems, one can formulate a complete qualitative dynamical
theory which includes notions of invariance, attractors, repellers, as well as connecting
orbits. The global dynamical behavior of a combinatorial multivector field is encoded in
the notion of a Morse decomposition, and it can be studied further using algebraic and
topological tools such as the Conley index [@conley:78a; @mischaikow:mrozek:02a]
and connection matrices [@franzosa:89a; @mrozek:wanner:25a]. Moreover, if the
combinatorial multivector field is generated from a classical dynamical system,
it is possible to derive statements about the underlying dyamics of the original
system [@mrozek:etal:22a; @thorpe:wanner:p24a; @thorpe:wanner:p24b]. In many cases,
these statements can eventually be turned into computer-assisted proofs
[@stephens:wanner:14a]. The Julia [@bezanson:17a] package `ConleyDynamics.jl`
[@conleydynamics] provides computational tools for combinatorial topological
dynamics, and should be of interest to both researchers and students which
are curious about this emerging field.

# Statement of need

To the best of our knowledge, `ConleyDynamics.jl` is the first public domain
software that provides a full implementation of all aspects of combinatorial
topological dynamics within a unified framework. Specifically, the following
three features form the core of the software:

* General Lefschetz complexes are implemented as the underlying combinatorial
  structure. They can be formulated over the field of rational numbers, or over
  any finite field of prime order. In addition, `ConleyDynamics.jl` provides
  specialized functions for simplicial and cubical complexes, which are two
  important subtypes of Lefschetz complexes.
* To describe the dynamics, general multivector fields are supported. As a
  special case these include Forman combinatorial vector fields. Functions
  are provided to generate such multivector fields from planar and
  three-dimensional autonomous ordinary differential equations. In addition,
  one can analyze the global dynamics of a combinatorial dynamical system using
  functions that determine Morse decompositions and Morse intervals.
* In order to implement all major aspects of Conley theory, algebraic
  topological tools are needed as well. Since the Conley index is given by
  homology groups, the package provides functions for computing regular and
  relative homology of Lefschetz complexes. In addition, it is possible to
  determine the connection matrix associated with a Morse decomposition.
  All of these computations can be performed over any of the fields
  specified above.

In the following we briefly describe each of these three features in more
detail.

General Lefschetz complexes as defined in [@lefschetz:42a] have usually
been passed over in favor of their two celebrated special subtypes: Simplicial
and cubical complexes. Simplicial complexes [@edelsbrunner:harer:10a;
@munkres:84a] form the foundation for many computational tools, such as
for example triangular meshes. In some cases, cubical complexes
[@kaczynski:etal:04a] are more convenient, such as in the analysis
of images or voxel data. However, for the analysis of classical dynamics
both of these combinatorial structures have disadvantages, as indicated
for example in [@boczko:etal:07a]. While one suitable generalization could
be CW complexes [@massey:91a; @dlotko:etal:11a], the package `ConleyDynamics.jl`
provides the most general structure, a Lefschetz complex. To keep the software
framework self-contained, special functions for the generation of both simplicial
and cubical complexes have been implemented as well. Note that even though parts
of the theory were formulated in the context of general finite topological
spaces [@alexandrov:37a; @lipinski:etal:23a], these are not provided by
`ConleyDynamics.jl`. The main reason for this is the lack of a notion
of connection matrix in this situation.

We now turn to the description of dynamics on combinatorial structures.
In the years since Forman introduced combinatorial vector fields, they
have found numerous applications in areas such as visualization and mesh
compression, graph braid groups, homology computation, astronomy, the study
of Cech and Delaunay complexes, and many others. See [@batko:etal:20a] for
detailed references. Despite these applications, there is no general purpose
software for working with such combinatorial vector fields. While Forman
vector fields have close ties with classical dynamics [@kaczynski:etal:16a;
@batko:etal:20a; @mrozek:wanner:21a], they are of limited use for the analysis
of a given classical flow. For example, in a combinatorial Forman vector field
chaos can only be generated through the choice of an appropriate underlying
combinatorial structure [@mrozek:etal:22a]. For this reason, the notion of
combinatorial multivector fields was introduced in [@mrozek:17a;
@lipinski:etal:23a], and it allows for the full spectrum of dynamical
behavior.

![A sample combinatorial multivector field on a small simplicial complex, together
with its Morse decomposition and Morse set Conley indices.\label{mvffig}](multivectorex.png){ width=90% }

A sample multivector field is shown in \autoref{mvffig}. Losely speaking, a
combinatorial multivector field partitions the cells of the underlying Lefschetz
complex into so-called locally closed subsets. (Precise definitions can be found
in [@conleydynamics; @lipinski:etal:23a].) Within each multivector, any
flow direction is admissible, as well as any flow towards the boundary of a
cell. In the figure, this indicates that there are stationary states, or 
equilibrium solutions, which give rise to constant invariant dynamics. In 
addition, one can observe periodic and more general recurrent behavior, as
well as transitions between invariant sets. The overall dynamics is described
via a Morse decomposition, which identifies essential invariant sets, gives
their homological Conley index, and describes transitional behavior between
these Morse sets, see for example the right-most panel in \autoref{mvffig}.

The package `ConleyDynamics.jl` provides functions which produce the Morse
decomposition for any given combinatorial multivector field, on any underlying
Lefschetz complex. These computations are fast, as they rely on the Julia
package `Graphs.jl` [@Graphs2021]. In addition, `ConleyDynamics.jl` allows for
the automatic construction of combinatorial multivector fields based on the 
concept of dynamical transitions [@mrozek:wanner:25a; @thorpe:wanner:p24a;
@thorpe:wanner:p24b]. This amounts to only having to specify essential dynamical
behavior, and then constructing the unique smallest multivector field which
realizes this behavior. The construction again relies on `Graphs.jl`,
and has been parallelized. Finally, functions are provided that use this
approach to create a combinatorial multivector field from a classical planar
or three-dimensional ordinary differential equation, on any underlying
Lefschetz complex that discretizes the relevant part of phase space. This
leads to a considerably more flexible discretization than the one described
in [@boczko:etal:07a]. Notice also that there is no publicly available
implementation of their algorithm.

We now address the algebraic topological aspects of `ConleyDynamics.jl`.
For the homology computations over the field types specified earlier,
we implemented the classical algorithm described in [@edelsbrunner:harer:10a].
There is better-performing software available for the computation of homology,
such as for example [@gudhi:24a]. However, our goal was the creation of a
unified framework for combinatorial topological dynamics, which allows for
general fields and Lefschetz complexes. It was therefore easier to include
the classical algorithm, which is more than sufficient for the purposes
of `ConleyDynamics.jl`.

The situation is different for the second algebraic topological aspect,
the computation of connection matrices. The first algorithm for the
computation of connection matrices was described in [@harker:etal:21a].
Their software, however, is geared towards their specific applications.
For this reason, `ConleyDynamics.jl` implements the recent method
described in [@dey:etal:24a], which leads to a persistence-like 
algorithm for the computation of connection matrices. Our implementation
extends their method to general fields and Lefschetz complexes. To the
best of our knowledge, no other open source implementations are available
at the present time.

# Examples

The usage of `ConleyDynamics.jl` is described in the extensive documentation
associated with [@conleydynamics]. We therefore only include a very brief
example and a few sample computations which analyze ordinary
differential equations. For the example shown in \autoref{mvffig},
the complete analysis of its global dynamics can be accomplished
with the following commands:

```julia
using ConleyDynamics
labels = ["A","B","C","D"]
simplices = [["A","B","C"],["B","C","D"]]
sclogo = create_simplicial_complex(labels,simplices)
mvflogo = [["A","AB"],["C","AC"],["B","BC","BD","BCD"]]
mdlogo = morse_sets(sclogo, mvflogo)
cmlogo = connection_matrix(sclogo, mvflogo)
```

The first command loads the package, and the next three initialize the
underlying Lefschetz complex. After specifying the multivector field in
the next line, the last two commands determine the Morse decomposition
and connection matrix, respectively. For more details, we refer the
reader to [@conleydynamics].

As mentioned earlier, `ConleyDynamics.jl` can determine a combinatorial
multivector field which encapsulates the possible dynamical behavior of
a classical autonomous ordinary differential equation. All that is needed
is the vector field describing the flow, and a Lefschetz complex discretization
of the relevant portion of phase space. In \autoref{flowfig1} we show the Morse
decompositions for two sample planar systems. They are based on underlying
phase space discretizations which are given by random Delaunay 
triangulations. The left system is a gradient system, which contains nine
equilibrium solutions. Enclosures for these are shown in \autoref{flowfig1},
and one can also easily determine their Conley indices, which describe their
stability properties in more detail. The right panel of \autoref{flowfig1}
depicts the Morse decomposition of a planar system with two periodic solutions,
which circle around an unstable equilibrium in their center. Also in this 
case, `ConleyDynamics.jl` is able to provide a description of the global
dynamics, including Conley indices and the connection matrix.

![Morse decompositions for two planar flows. In both cases, the underlying
Lefschetz complex is a random Delaunay triangulation. The identified Morse 
sets are shown in different colors.\label{flowfig1}](flowfig1.png)

![Morse intervals for a planar gradient flow. The yellow regions provide
encloures for heteroclinic orbits in the system.\label{flowfig2}](flowfig2.png)

Furthermore, `ConleyDynamics.jl` can provide additional information on the global
dynamics of the system. For the gradient differential equation underlying
the left image of \autoref{flowfig1}, it is possible to determine enclosures
for the connecting orbits between stationary states. These are shown in
\autoref{flowfig2}. In the left image, enclosures for the heteroclinic
orbits between four unstable equilibria and the four asymptotically stable
ones are depicted, while the image on the right illustrates the location
of four connecting orbits between the unstable solution at the center
and the other four unstable stationary states. For more details, see
[@conleydynamics].

As our final example we briefly illustrate the analysis of a three-dimensional
autonomous ordinary differential equation. In the left panel of \autoref{flowfig3}
we show the Morse decomposition of a spatial flow with seven equilibrium
solutions. The depicted enclosures are based on a coarse cubical grid, yet 
nevertheless they correctly identify the regions of interest. It is also
possible to determine stability information via the Conley index, and the panel
on the right illustrates connecting orbits between an unstable equilibrium
and two asymptotically stable ones.

![Morse decomposition and sample Morse interval for a three-dimensional
ordinary differential equation. The underlying Lefschetz complex is a
coarse cubical grid.\label{flowfig3}](flowfig3.png)

The above examples provide a brief glimpse into some of the capabilities of
`ConleyDynamics.jl`. For more details on the equations underlying these
examples, as well as numerous additional demonstrations and a full description
of all implemented functionality, we refer the reader to the manual which
accompanies [@conleydynamics].

# Acknowledgements

T.W. was partially supported by the Simons Foundation under Award 581334.

# References

