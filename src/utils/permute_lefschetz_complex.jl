export permute_lefschetz_complex

"""
    permute_lefschetz_complex(lc::LefschetzComplex,
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

    ncells1     = lc.ncells
    lcdim1      = lc.dim
    boundary1   = lc.boundary
    labels1     = lc.labels
    dimensions1 = lc.dimensions

    # Create the permuted data

    ncells2     = ncells1
    lcdim2      = lcdim1
    labels2     = labels1[perm]
    dimensions2 = dimensions1[perm]
    boundary2   = sparse_permute(boundary1, perm, perm)

    # Create the permuted dictionary
    
    indices2 = Dict{String,Int}([(labels2[k],k) for k in 1:length(labels2)])

    # Create the permuted Lefschetz complex

    lc2 = LefschetzComplex(ncells2, lcdim2, boundary2,
                           labels2, indices2, dimensions2)
    
    # Return the permuted Lefschetz complex

    return lc2
end

