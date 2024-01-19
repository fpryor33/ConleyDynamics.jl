export connection_matrix

"""
    cm = connection_matrix(lc, mvf; p=2)
    cm, cmbasis = connection_matrix(lc, mvf; p=2, returnbasis=true)

Compute a connection matrix for the multivector field `mvf` on the
Lefschetz complex `lc` over a finite field with `p` elements.

The arguments are typed as `lc::LefschetzComplex` and `mvf::Vector{Vector{Int}}`,
and the return object is of type `ConleyMorseCM`. If the optional argument
`returnbasis::Bool=true` is given, then the function returns a dictionary 
which gives the basis for the connection matrix columns in terms of the
original labels. If `p` is omitted, then `p=2` is used.
"""
function connection_matrix(lc::LefschetzComplex, mvf::Vector{Vector{Int}};
                           p::Int=2, returnbasis::Bool=false)
    #
    # Compute the connection matrix
    #

    # Copy boundary and label data
    
    bndmatrix = lc.boundary
    labels    = lc.label
    ppoly     = lc.poincare

    # Find an admissible order

    adorder, adbnd, psetvec, scc = admissible_order(bndmatrix, mvf)

    # Convert the boundary matrix to finite field format
    # For now we hardcode the characteristic of the finite field

    bndA = convert_matrix_gfp(adbnd,p)

    # Compute the connection matrix in reordered form

    cmRedP = deepcopy(bndA)
    if returnbasis
        cmMatrP, cmColsP, cmBasisP = cm_create!(cmRedP, psetvec; returnbasis=true)
    else
        cmMatrP, cmColsP = cm_create!(cmRedP, psetvec)
    end

    # Compute first few return variables
    
    cmMatr      = cmMatrP
    cmCols      = adorder[cmColsP]
    cmPoset     = psetvec[cmColsP]
    cmLabels    = labels[cmCols]

    # If desired, determine the basis dictionary

    if returnbasis
        basisdict = Dict{String,Vector{String}}()
        for k=1:length(cmColsP)
            basisStr = labels[cmCols[k]]
            basisVec = Vector{String}()
            for m=size(bndmatrix,2):-1:1
                if !(cmBasisP[m,cmColsP[k]]==0)
                    push!(basisVec,labels[adorder[m]])
                end
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

    if returnbasis
        return cm, basisdict
    else
        return cm
    end
end

"""
    cm = connection_matrix(lc, mvf; p=2)
    cm, cmbasis = connection_matrix(lc, mvf; p=2, returnbasis=true)

Compute a connection matrix for the multivector field `mvf` on the
Lefschetz complex `lc` over a finite field with `p` elements.

The arguments are typed as `lc::LefschetzComplex` and `mvf::Vector{Vector{String}}`,
and the return object is of type `ConleyMorseCM`. If the optional argument
`returnbasis::Bool=true` is given, then the function returns a dictionary 
which gives the basis for the connection matrix columns in terms of the
original labels. If `p` is omitted, then `p=2` is used.
"""
function connection_matrix(lc::LefschetzComplex, mvf::Vector{Vector{String}};
                           p::Int=2, returnbasis::Bool=false)
    #
    # Compute the connection matrix
    #
    # This is an alternative method for multivector fields with
    # a String representation.
    #
    newmvf = convert_mvf(mvf, lc)
    if returnbasis
        cm, cmbasis = connection_matrix(lc, newmvf; p=p, returnbasis=true)
        return cm, cmbasis
    else
        cm = connection_matrix(lc, newmvf; p=p)
        return cm
    end
end

