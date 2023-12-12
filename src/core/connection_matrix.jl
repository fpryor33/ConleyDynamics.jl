export connection_matrix

"""
    cm = connection_matrix(lc, mvf)

Compute a connection matrix for the multivector field `mvf` on the
Lefschetz complex `lc`.

The arguments are typed as `lc::LefschetzComplex` and `mvf::Vector{Vector{Int}}`,
and the return object is of type `ConleyMorseCM`.
"""
function connection_matrix(lc::LefschetzComplex, mvf::Vector{Vector{Int}})
    #
    # Compute the connection matrix
    #

    # Copy boundary and label data
    
    bndmatrix = lc.boundary
    labels    = lc.label
    ppoly     = lc.poincare

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

    # Determine the Morse sets and their Poincare polynomials
    # using the original labels

    poset_indices = sort(union(cmPoset))
    cmMorseSets = Vector{Vector{String}}()
    for k = 1:length(poset_indices)
        push!(cmMorseSets,labels[scc[poset_indices[k]]])
    end

    cmPoincare  = Vector{typeof(ppoly[1])}()
    for k=1:length(poset_indices)
        ppolytemp = 0 * ppoly[1]
        for j=1:length(cmPoset)
            if cmPoset[j] == poset_indices[k]
                ppolytemp = ppolytemp + ppoly[cmCols[j]]
            end
        end
        push!(cmPoincare,ppolytemp)
    end

    # Renumber the poset indices in consecutive order

    renumber_poset!(cmPoset)

    # Return the connection matrix information
    
    cm = ConleyMorseCM{typeof(cmMatr),typeof(ppoly[1])}(
                cmMatr, cmCols, cmPoset, cmLabels, cmMorseSets, cmPoincare)

    return cm
end

