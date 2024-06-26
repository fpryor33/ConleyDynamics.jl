export connection_matrix

"""
    connection_matrix(lc::LefschetzComplex, mvf::MultiVectorField;
                      [returnbasis::Bool])

Compute a connection matrix for the multivector field `mvf` on the
Lefschetz complex `lc` over the field associated with the Lefschetz
complex boundary matrix.

The function returns an object of type `ConleyMorseCM`. If the optional
argument `returnbasis::Bool=true` is given, then the function also returns
a dictionary which gives the basis for the connection matrix columns in
terms of the original labels.
"""
function connection_matrix(lc::LefschetzComplex, mvfarg::MultiVectorField;
                           returnbasis::Bool=false)
    #
    # Compute the connection matrix
    #

    # Convert the multivector field to Int if necessary

    if mvfarg isa Vector{Vector{String}}
        mvf = convert_mvf(lc, mvfarg)
    else
        mvf = mvfarg
    end

    # Copy boundary and label data
    
    lcdim     = lc.dim
    bndmatrix = lc.boundary
    labels    = lc.labels
    celldims  = lc.dimensions

    # Find an admissible order

    adorder, adbnd, psetvec, scc = admissible_order(bndmatrix, mvf)

    # Compute the connection matrix in reordered form

    cmRedP = deepcopy(adbnd)
    if returnbasis
        cmMatrP, cmColsP, cmBasisP = cm_reduce!(cmRedP, psetvec; returnbasis=true)
    else
        cmMatrP, cmColsP = cm_reduce!(cmRedP, psetvec)
    end

    # Compute first few return variables
    
    cmMatr      = cmMatrP
    cmCols      = adorder[cmColsP]
    cmPoset     = psetvec[cmColsP]
    cmLabels    = labels[cmCols]

    # If desired, determine the basis dictionary

    if returnbasis
        basisdict = Dict{String,Vector{String}}()
        for k=1:length(cmBasisP)
            basisStr = labels[cmCols[k]]
            basisVec = Vector{String}()
            for m=1:length(cmBasisP[k])
                push!(basisVec,labels[adorder[cmBasisP[k][m]]])
            end
            basisdict[basisStr] = basisVec
        end
    end

    # Determine the Morse sets and their Poincare polynomials
    # using the original labels

    poset_indices = sort(union(cmPoset))
    cmMorseSets = Vector{Vector{String}}()
    for k = 1:length(poset_indices)
        push!(cmMorseSets,labels[scc[poset_indices[k]]])
    end

    cmPoincare  = Vector{Vector{Int}}()
    for k=1:length(poset_indices)
        ppolytemp = fill(Int(0),lcdim+1)
        for j=1:length(cmPoset)
            if cmPoset[j] == poset_indices[k]
                ppolytemp[celldims[cmCols[j]]+1] += 1
            end
        end
        push!(cmPoincare,ppolytemp)
    end

    # Renumber the poset indices in consecutive order

    renumber_poset!(cmPoset)

    # Construct the Conley complex as a Lefschetz complex

    cc_ncells  = length(cmLabels)
    cc_bnd     = cmMatr
    cc_labels  = cmLabels
    cc_dims    = Vector{Int}([lc.dimensions[lc.indices[cmLabels[k]]]
                              for k in 1:cc_ncells])
    cc_dim     = maximum(cc_dims)
    cc_indices = Dict{String,Int}(cc_labels[j] => j
                                   for j=1:length(cc_labels))

    # Create the new Lefschetz complex and return it

    cmCC = LefschetzComplex(cc_ncells, cc_dim, cc_bnd,
                            cc_labels, cc_indices, cc_dims)
 
    # Return the connection matrix information
    
    cm = ConleyMorseCM{typeof(cmMatr.zero)}(
                cmMatr, cmCols, cmPoset, cmLabels, cmMorseSets, cmPoincare, cmCC)

    if returnbasis
        return cm, basisdict
    else
        return cm
    end
end

