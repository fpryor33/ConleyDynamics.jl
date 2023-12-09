export connection_matrix

"""
    cmMatr, cmCols, cmPoset, cmMorseSets, cmLabels = connection_matrix(lc, mvf)

Compute a connection matrix for the multivector field `mvf` on the
Lefschetz complex `lc`.

The arguments are typed as `lc::LefschetzComplex` and `mvf::Vector{Vector{Int}}`.

# Return Arguments
* `cmMatr`: connection matrix
* `cmCols`: columns from `bndmatrix` associated with columns of `cmMatr`
* `cmPoset`: poset labels for the Morse sets
* `cmMorseSets`: Morse sets
* `cmLabels`: labels from `labels` corresponding to columns of `cmMatr`
"""
function connection_matrix(lc::LefschetzComplex, mvf::Vector{Vector{Int}})
    #
    # Compute the connection matrix
    #

    # Copy boundary and label data
    
    bndmatrix = lc.boundary
    labels    = lc.label

    # Find an admissible order

    adorder, psetvec, scc = admissible_order(bndmatrix, mvf)

    # Convert the boundary matrix to finite field format
    # For now we hardcode the characteristic of the finite field

    p = 2
    bndA = convert_matrix_gfp(bndmatrix[adorder,adorder],p)

    # Compute the connection matrix in reordered form

    cmMatrP, cmColsP, cmRedP = cm_create(bndA, psetvec)

    # Compute first few return variables
    
    cmMatr      = cmMatrP
    cmCols      = adorder[cmColsP]
    cmPoset     = psetvec[cmColsP]
    cmLabels    = labels[cmCols]

    # Determine the Morse sets using the original labels

    poset_indices = sort(union(cmPoset))
    cmMorseSets = Vector{Vector{typeof(labels[1])}}()
    for k = 1:length(poset_indices)
        push!(cmMorseSets,labels[scc[poset_indices[k]]])
    end

    # Renumber the poset indices in consecutive order

    renumber_poset!(cmPoset)

    # Return the connection matrix information

    return cmMatr, cmCols, cmPoset, cmMorseSets, cmLabels
end

