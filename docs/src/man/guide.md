# Tutorial

In this tutorial we explain the basic usage of the main components
of ConleyDynamics. It is not meant to be exhaustive, since more 
details will be provided in the more indiviualized sections. Also,
precise mathematical definitions will be delayed until then.

## Creating Lefschetz Complexes

The fundamental mathematical object for ConleyDynamics is a Lefschetz
complex. For now we note that both simplicial complexes and cubical
complexes are special cases, and ConleyDynamics provides convenient
interfaces for generating them.

We begin by considering the case of a simplicial complex. Recall that
an *abstract simplicial complex* ``K`` is just a collection of finite
sets, called *simplices*, which is closed under taking subsets. In
other words, every subset of a simplex is again a simplex. Each
simplex has an associated *dimension*, which is one less than the 
number of its elements. 

```@example T1
using ..ConleyDynamics # hide
labels = ["A","B","C","D","E","F"]
simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
sc = create_simplicial_complex(labels,simplices)
fieldnames(typeof(sc))
```


```@example T1
println(sc.ncells)
```


```@example T1
println(sc.dim)
```


```@example T1
println(sc.boundary)
```


```@example T1
println(sc.labels)
```

```@example T1
println(sc.dimensions)
```



```@example T1
println(sc.indices)
```







## Computing Homology and Persistence



```@example T1
homology(sc,p=0)
```

```@example T1
relative_homology(sc, [1,6], p=0)
```

```@example T1
filtration = [1,1,1,2,2,2,1,1,1,3,2,2,2,4]
phsingles, phpairs = persistent_homology(sc, filtration, p=2)
```



## Finding Connection Matrices


![The logo multivector field](img/multivectorex.png)

We begin by discussion an elementary example from Mrozek & Wanner[^1].

[^1]:
    Marian Mrozek, Thomas Wanner: [Connection matrices in combinatorial
    topological dynamics](https://arxiv.org/abs/2103.04269),
    *Preprint*, submitted for publication, 115 pp, 2023.

## Working with Sparse Matrices



