export simplicial_torus
export simplicial_klein_bottle
export simplicial_projective_plane

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
    # Create a triangulation of the projective_plane
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

