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

---

# Summary

Combinatorial topological dynamics is concerned with the qualitative study of
dynamical behavior on discrete combinatorial structures. It was originally developed
in the context of combinatorial vector fields [@forman:98a; @forman:98b], and has since
been extended to combinatorial multivector fields [@mrozek:17a; @lipinski:etal:23a] on
Lefschetz complexes [@lefschetz:42a]. For such systems, one can formulate a complete
qualitative theory which includes notions of invariance, attractors, repellers, and
connecting orbits. The global dynamical behavior is encoded in a Morse decomposition,
and it can be studied further using algebraic topological tools such as the Conley
index [@conley:78a; @mischaikow:mrozek:02a] and connection matrices [@franzosa:89a;
@mrozek:wanner:25a]. If the combinatorial multivector field is generated from a
classical flow, one can derive statements about the underlying dyamics of the original
system [@mrozek:etal:22a; @thorpe:wanner:p24a; @thorpe:wanner:p24b]. The Julia
[@bezanson:17a] package `ConleyDynamics.jl` provides computational tools for
combinatorial topological dynamics, and should be of interest to both
researchers and students which are curious about this emerging field.

# Statement of need

`ConleyDynamics.jl` provides a full implementation of all aspects of
combinatorial topological dynamics within a unified framework.
The following form the core of the software:

* General Lefschetz complexes are implemented as the underlying combinatorial
  structure. They can be formulated over the field of rational numbers, or over
  any finite field of prime order. Specialized functions for simplicial and
  cubical complexes are provided.
* To describe the dynamics, general multivector fields are supported, which
  include Forman combinatorial vector fields. Functions are provided to generate
  multivector fields from planar and three-dimensional autonomous ordinary
  differential equations. One can analyze the global dynamics of a combinatorial
  dynamical system by determining Morse decompositions and Morse intervals.
* In order to implement all major aspects of Conley theory, algebraic
  topological tools are needed as well. The package provides functions
  for computing regular and relative homology of Lefschetz complexes, as
  well as the connection matrix associated with a Morse decomposition.

General Lefschetz complexes [@lefschetz:42a] have usually
been passed over in favor of their two celebrated special subtypes: Simplicial
and cubical complexes. Simplicial complexes [@edelsbrunner:harer:10a;
@munkres:84a] form the foundation for many computational tools, such as
for example triangular meshes. Sometimes cubical complexes
[@kaczynski:etal:04a] are more convenient, such as in the analysis
of images or voxel data. However, for the analysis of classical dynamics
both of these combinatorial structures have disadvantages, as indicated
in [@boczko:etal:07a]. While one suitable generalization could be CW complexes
[@massey:91a; @dlotko:etal:11a], `ConleyDynamics.jl` provides the most general
structure, a Lefschetz complex. To keep the software framework self-contained,
special functions for simplicial and cubical complexes have been implemented
as well.

Since they were first introduced by Forman, combinatorial vector fields have
found numerous applications in areas such as visualization and mesh
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

`ConleyDynamics.jl` provides functions which analyze the global dynamics
of a combinatorial multivector field, on any underlying Lefschetz complex,
in terms of its Morse decomposition. These computations are fast, as they
rely on `Graphs.jl` [@Graphs2021]. Combinatorial multivector fields can be
constructed based on the concept of dynamical transitions [@mrozek:wanner:25a;
@thorpe:wanner:p24a; @thorpe:wanner:p24b]. Finally, functions are provided
that use this approach to create a combinatorial multivector field from a
classical planar or three-dimensional ordinary differential equation, on
any underlying Lefschetz complex that discretizes the relevant part of
phase space. This leads to more flexible discretizations than the ones
described in [@boczko:etal:07a].

For the necessary homology computations, the classical algorithm described
in [@edelsbrunner:harer:10a] is used. Its implementation in
`ConleyDynamics.jl` allows for computations over general prime finite fields
or over the rationals, and one can work with arbitrary Lefschetz complexes. 
Note, however, that there is better-performing software available, such as
for example [@gudhi:24a]. 

The first algorithm for the computation of connection matrices was
described in [@harker:etal:21a]. Their software, however, is geared
towards specific applications. For this reason, `ConleyDynamics.jl`
implements the recent method described in [@dey:etal:24a]. Our
implementation extends their method to general fields and
Lefschetz complexes.

![Two sample planar flows. The left image shows a gradient system with
nine equilibrium solutions and many connecting orbits betweeen them. The 
system depicted on the right has one stable equilibrium at the origin and
two surrounding periodic orbits.\label{flowfig0}](flowfig0.png)

# Examples

`ConleyDynamics.jl` can determine a combinatorial multivector field which
encapsulates the possible dynamical behavior of a classical flow, based on
the underlying vector field and a Lefschetz complex discretization of the
relevant portion of phase space. 

![Morse decompositions for the two planar flows shown in the previous figure.
In both cases, the underlying Lefschetz complex is a random Delaunay triangulation.
The identified Morse sets are shown in different colors.\label{flowfig1}](flowfig1.png)

![Morse intervals for the planar gradient flow. The yellow regions provide
enclosures for heteroclinic orbits in the system.\label{flowfig2}](flowfig2.png)

Consider for example the two planar flows
depicted in \autoref{flowfig0}. For these systems, \autoref{flowfig1} shows
Morse decompositions based on underlying random Delaunay triangulations.
The left system contains nine equilibrium solutions, within the colored
enclosures. The right panel of \autoref{flowfig1} depicts the Morse decomposition
of a planar system with two periodic solutions, which circle around an
equilibrium in their center. In both cases, `ConleyDynamics.jl` also
provides Conley indices and the connection matrix. In addition, one can
determine additional information on the global dynamics of these system.
For example, \autoref{flowfig2} visualizes enclosures for the connecting
orbits between stationary states, for the system in the left panel of
\autoref{flowfig1}. For more details and examples, see the extensive
documentation accompanying [@conleydynamics], and the two recent papers
[@thorpe:wanner:p24a] and [@thorpe:wanner:p24b]. In addition, all examples 
in the book [@mrozek:wanner:25a] were computed using `ConleyDynamics.jl`.

# Acknowledgements

T.W. was partially supported by the Simons Foundation under Award 581334.

# References

