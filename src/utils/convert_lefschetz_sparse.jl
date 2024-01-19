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

    ncells1   = lc.ncells
    boundary1 = lc.boundary
    label1    = lc.label
    index1    = lc.index
    poincare1 = lc.poincare

    # Create the sparse boundary matrix

    boundary2 = sparse_from_full(boundary1)

    # Create and return the sparse Lefschetz complex

    lc2 = LefschetzComplex{typeof(poincare1[1])}(ncells1,
                           boundary2, label1, index1, poincare1)
    return lc2
end

