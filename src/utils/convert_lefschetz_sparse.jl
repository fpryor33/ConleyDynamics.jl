export convert_lefschetz_sparse

"""
    lc2 = convert_lefschetz_sparse(lc::LefschetzComplex)

Convert the boundary matrix of a Lefschetz complex to sparse.
"""
function convert_lefschetz_sparse(lc::LefschetzComplex)
    #
    # Convert a Lefschetz complex to one with sparse boundary
    #

    # Copy the data from the original complex

    ncells1     = lc.ncells
    dim1        = lc.dim
    boundary1   = lc.boundary
    labels1     = lc.labels
    indices1    = lc.indices
    dimensions1 = lc.dimensions

    # Create the sparse boundary matrix

    boundary2 = sparse_from_full(boundary1)

    # Create and return the sparse Lefschetz complex

    lc2 = LefschetzComplex(ncells1, dim1, boundary2,
                           labels1, indices1, dimensions1)
    return lc2
end

