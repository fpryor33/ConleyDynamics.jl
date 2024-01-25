export example_MW_fig04

"""
    lc1, lc2, mvf = example_MW_fig04()

Create two representations of the Lefschetz complex and
the multivector field for the example from Figure 4 in
the connection matrix paper by *Mrozek & Wanner*.

The complexes `lc1` and `lc2` are just two representations
of the same complex, but they lead to different connection
matrices.

# Examples
```jldoctest
julia> lc1, lc2, mvf = example_MW_fig04();

julia> cm1 = connection_matrix(lc1, mvf);

julia> cm2 = connection_matrix(lc2, mvf);

julia> cm1.cm
4×4 Matrix{Nemo.fpFieldElem}:
 0  0  0  0
 0  0  0  1
 0  0  0  1
 0  0  0  0

julia> cm2.cm
4×4 Matrix{Nemo.fpFieldElem}:
 0  0  0  0
 0  0  0  0
 0  0  0  1
 0  0  0  0
```
"""
function example_MW_fig04()
    #
    # Create two representations of the Lefschetz complex
    # and the multivector field for the example from 
    # Figure 4 in the connection matrix paper
    # by *Mrozek & Wanner*.
    
    nc = 6

    # Create the vector of labels

    labelvec = Vector{String}()
    push!(labelvec,"A","B","a","b","c","alpha")

    # Create the label to index dictionary

    indexdict = Dict{String,Int}([(labelvec[k],k) for k in 1:length(labelvec)])

    # Create the vector of simplex dimensions
    
    sdvec = [0, 0, 1, 1, 1, 2]

    # Create the boundary matrix

    bndmatrix = zeros(Int, nc, nc)

    # Start with the edges

    bndmatrix[[1,2],3] = [1; 1]     # a
    bndmatrix[[1,2],4] = [1; 1]     # b
    bndmatrix[[1,2],5] = [1; 1]     # c

    # Move on to the 2-cell

    bndmatrix[[3,4],6] = [1; 1]     # alpha 

    # Construct the Lefschetz complex struct
    
    lc = LefschetzComplex(nc, Int(2), bndmatrix,
                          labelvec, indexdict, sdvec)

    # Create a second version of the Lefschetz complex via permutation

    perm = [1, 2, 5, 3, 4, 6]      # Change a, b, c to the order c, a, b
    labelvec2  = labelvec[perm]
    sdvec2     = sdvec[perm]
    bndmatrix2 = bndmatrix[perm,perm]
    indexdict2 = Dict{String,Int}([(labelvec2[k],k) for k in 1:length(labelvec2)])

    lc2 = LefschetzComplex(nc, Int(2), bndmatrix2,
                           labelvec2, indexdict2, sdvec2)

    # Create the common part of the combinatorial vector fields
    
    mvf = Vector{Vector{String}}()
    push!(mvf,["A", "a"])       # A - a
    push!(mvf,["B", "c"])       # B - c

    # Return the example data

    return lc, lc2, mvf
end

