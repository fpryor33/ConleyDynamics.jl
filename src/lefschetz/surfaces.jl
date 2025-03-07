export simplicial_torus
export simplicial_klein_bottle
export simplicial_projective_plane
export simplicial_torsion_space

"""
    sc = simplicial_torus(p::Int)

Create a triangulation of the two-dimensional torus.

The function returns a simplicial complex which represents a
two-dimensional torus. The argument `p` specifies the characteristic
of the underlying field. This triangulation is taken from Figure 6.4
in Munkres' book on Algebraic Topology. The boundary vertices are
labeled as letters as in the book, the five center vertices are labeled
by `1` through `5`.

# Examples
```jldoctest
julia> println(homology(simplicial_torus(0)))
[1, 2, 1]

julia> println(homology(simplicial_torus(2)))
[1, 2, 1]

julia> println(homology(simplicial_torus(3)))
[1, 2, 1]
```
"""
function simplicial_torus(p::Int)
    # 
    # Create a triangulation of the torus
    #
    
    # Create the simplicial complex

    labels = ["a","b","c","d","e","1","2","3","4","5"]
    simplices = [["1","2","5"],["2","3","5"],["3","4","5"],["1","4","5"],
                 ["a","b","1"],["a","d","1"],["b","c","1"],["c","1","2"],
                 ["a","c","2"],["a","d","2"],["d","2","3"],["d","e","3"],
                 ["e","a","3"],["a","c","3"],["b","c","3"],["b","3","4"],
                 ["a","e","4"],["a","b","4"],["d","e","1"],["e","1","4"]]
    lc = create_simplicial_complex(labels, simplices, p=p)
    return lc
end

"""
    sc = simplicial_klein_bottle(p::Int)

Create a triangulation of the two-dimensional Klein bottle.

The function returns a simplicial complex which represents the
two-dimensional Klein bottle. The argument `p` specifies the
characteristic of the underlying field. This triangulation is
taken from Figure 6.6 in Munkres' book on Algebraic Topology.
The boundary vertices are labeled as letters as in the book,
the five center vertices are labeled by `1` through `5`.

# Examples
```jldoctest
julia> println(homology(simplicial_klein_bottle(0)))
[1, 1, 0]

julia> println(homology(simplicial_klein_bottle(2)))
[1, 2, 1]

julia> println(homology(simplicial_klein_bottle(3)))
[1, 1, 0]
```
"""
function simplicial_klein_bottle(p::Int)
    # 
    # Create a triangulation of the Klein bottle
    #
    
    # Create the simplicial complex

    labels = ["a","b","c","d","e","1","2","3","4","5"]
    simplices = [["1","2","5"],["2","3","5"],["3","4","5"],["1","4","5"],
                 ["a","b","1"],["a","d","1"],["b","c","1"],["c","1","2"],
                 ["a","c","2"],["a","e","2"],["e","2","3"],["d","e","3"],
                 ["d","a","3"],["a","c","3"],["b","c","3"],["b","3","4"],
                 ["a","e","4"],["a","b","4"],["d","e","1"],["e","1","4"]]
    lc = create_simplicial_complex(labels, simplices, p=p)
    return lc
end

"""
    sc = simplicial_projective_plane(p::Int)

Create a triangulation of the projective plane.

The function returns a simplicial complex which represents the
projective plane. The argument `p` specifies the characteristic
of the underlying field. This triangulation is taken from
Figure 6.6 in Munkres' book on Algebraic Topology. The boundary
vertices are labeled as letters as in the book, the five center
vertices are labeled by `1` through `5`.

# Examples
```jldoctest
julia> println(homology(simplicial_projective_plane(0)))
[1, 0, 0]

julia> println(homology(simplicial_projective_plane(2)))
[1, 1, 1]

julia> println(homology(simplicial_projective_plane(3)))
[1, 0, 0]
```
"""
function simplicial_projective_plane(p::Int)
    # 
    # Create a triangulation of the projective plane
    #
    
    # Create the simplicial complex

    labels = ["a","b","c","d","e","f","1","2","3","4","5"]
    simplices = [["1","2","5"],["2","3","5"],["3","4","5"],["1","4","5"],
                 ["a","b","1"],["a","f","1"],["b","c","1"],["c","1","2"],
                 ["c","d","2"],["d","e","2"],["e","2","3"],["e","f","3"],
                 ["f","a","3"],["a","b","3"],["b","c","3"],["c","3","4"],
                 ["c","d","4"],["d","e","4"],["e","f","1"],["e","1","4"]]
    lc = create_simplicial_complex(labels, simplices, p=p)
    return lc
end

"""
    sc = simplicial_torsion_space(n::Int, p::Int)

Create a triangulation of a space with 1-dimensional torsion.

The function returns a simplicial complex which has the following
integer homology groups:

- In dimension 0 it is the group of integers.
- In dimension 1 it is the integers modulo `n`.
- In dimension 2 it is the trivial group.

In other words, the simplicial complex has nontrivial torsion
in dimension 1. It is a triangulation of an `n`-gon, in which
all boundary edges are oriented counterclockwise, and all of
these edges are identified. The parameter `p` specifies the
characteristic of the underlying field.

# Examples
```jldoctest
julia> println(homology(simplicial_torsion_space(6,2)))
[1, 1, 1]

julia> println(homology(simplicial_torsion_space(6,3)))
[1, 1, 1]

julia> println(homology(simplicial_torsion_space(6,5)))
[1, 0, 0]
```
"""
function simplicial_torsion_space(n::Int, p::Int)
    # 
    # Create a triangulation of a space with 1-d torsion
    #
    
    # Make sure that n is at least 3

    if n < 3
        error("n has to be at least 3!")
    end

    # Compute the number of digits for the labels

    labwidth = Int(ceil(log(0.5+n)/log(10)))
    labformat = "%0" * string(labwidth) * "d%s"

    # Create the vector of labels

    labels = Vector{String}()

    vlabel = Printf.format(Printf.Format(labformat), 0, "c")   # Center point
    push!(labels, vlabel)
    vlabel = Printf.format(Printf.Format(labformat), 0, "w")   # Boundary point 1
    push!(labels, vlabel)
    vlabel = Printf.format(Printf.Format(labformat), 0, "x")   # Boundary point 2
    push!(labels, vlabel)
    vlabel = Printf.format(Printf.Format(labformat), 0, "y")   # Boundary point 3
    push!(labels, vlabel)
    vlabel = Printf.format(Printf.Format(labformat), 0, "z")   # Boundary point 4
    push!(labels, vlabel)

    for k = 1:n
        vlabel = Printf.format(Printf.Format(labformat), k, "s")
        push!(labels, vlabel)
        vlabel = Printf.format(Printf.Format(labformat), k, "t")
        push!(labels, vlabel)
    end

    # Create the vector of top-level triangles

    triangles = Vector{Vector{String}}()

    vlabelc = Printf.format(Printf.Format(labformat), 0, "c")
    for k = 0:n-1
        vlabels  = Printf.format(Printf.Format(labformat), k+1, "s")
        vlabelt  = Printf.format(Printf.Format(labformat), k+1, "t")
        kn = mod(k+1,n)
        vlabelsn = Printf.format(Printf.Format(labformat), kn+1, "s")
        push!(triangles, [vlabelc, vlabels, vlabelt])
        push!(triangles, [vlabelc, vlabelt, vlabelsn])
    end

    vlabelb1 = Printf.format(Printf.Format(labformat), 0, "w")
    vlabelb2 = Printf.format(Printf.Format(labformat), 0, "x")
    vlabelb3 = Printf.format(Printf.Format(labformat), 0, "y")
    vlabelb4 = Printf.format(Printf.Format(labformat), 0, "z")
    for k = 0:n-1
        vlabels  = Printf.format(Printf.Format(labformat), k+1, "s")
        vlabelt  = Printf.format(Printf.Format(labformat), k+1, "t")
        kn = mod(k+1,n)
        vlabelsn = Printf.format(Printf.Format(labformat), kn+1, "s")
        push!(triangles, [vlabels,  vlabelb1, vlabelb2])
        push!(triangles, [vlabels,  vlabelb2, vlabelb3])
        push!(triangles, [vlabels,  vlabelt,  vlabelb3])
        push!(triangles, [vlabelt,  vlabelb3, vlabelb4])
        push!(triangles, [vlabelt,  vlabelsn, vlabelb4])
        push!(triangles, [vlabelsn, vlabelb1, vlabelb4])
    end

    # Create the simplicial complex

    lc = create_simplicial_complex(labels, triangles, p=p)
    return lc
end

