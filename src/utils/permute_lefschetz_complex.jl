export permute_lefschetz_complex

"""
    lc2 = permute_lefschetz_complex(lc::LefschetzComplex,
                                    permutation::Vector{Int})

Permute the indices of a Lefschetz complex.

The vector `permutation` contains a permutation of the indices for the
given Lefschetz complex `lc`. If no permutation is specified, or if the
length of the vector is not correct, then a randomly generated one will
be used.
"""
function permute_lefschetz_complex(lc::LefschetzComplex,
                                   perm::Vector{Int}=Vector{Int}([]))
    #
    # Create a Lefschetz complex which permutes the cell indices
    #
    
    # Create a random permutation if none is set

    if !(length(perm) == lc.ncells)
        perm = randperm(lc.ncells)
    end

    # Copy the data from the original complex

    ncells1   = lc.ncells
    boundary1 = lc.boundary
    label1    = lc.label
    poincare1 = lc.poincare

    # Create the permuted data

    ncells2   = ncells1
    label2    = label1[perm]
    poincare2 = poincare1[perm]

    if typeof(boundary1) == Matrix{Int}
        boundary2 = boundary1[perm,perm]
    else
        boundary2 = sparse_permute(boundary1, perm, perm)
    end

    # Create the permuted dictionary
    
    index2 = Dict{String,Int}([(label2[k],k) for k in 1:length(label2)])

    # Create the permuted Lefschetz complex

    lc2 = LefschetzComplex{typeof(poincare2[1])}(ncells2,
                           boundary2, label2, index2, poincare2)
    
    # Return the permuted Lefschetz complex

    return lc2
end

