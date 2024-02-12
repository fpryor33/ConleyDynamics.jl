export connection_matrix

"""
    connection_matrix(lc::LefschetzComplex, mvf::MultiVectorField;
                      [p::Int,] [returnbasis::Bool])

Compute a connection matrix for the multivector field `mvf` on the
Lefschetz complex `lc` over a finite field with `p` elements, or 
over the rationals.

The function returns an object of type `ConleyMorseCM`. If the optional
argument `returnbasis::Bool=true` is given, then the function also returns
a dictionary which gives the basis for the connection matrix columns in
terms of the original labels. If `p` is omitted, then the Lefschetz
complex boundary has to has been specified over a field. If the boundary
matrix is an integer matrix, `p` has to be chosen.
"""
function connection_matrix(lc::LefschetzComplex, mvfarg::MultiVectorField;
                           p::Int=-1, returnbasis::Bool=false)
    #
    # Compute the connection matrix
    #

    # Convert the multivector field to Int if necessary

    if mvfarg isa Vector{Vector{String}}
        mvf = convert_mvf(mvfarg, lc)
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

    # Convert the boundary matrix to finite field or rational format

    if (adbnd isa SparseMatrix{Int}) && (adbnd.char == 0)
        if (p == -1)
            error("Homology over Z is not supported, specify p!")
        elseif (p >= 0)
            bndA = convert_matrix_gfp(adbnd, p)
        else
            error("Wrong characteristic p!")
        end
    else
        bndA = adbnd
        if !(adbnd.char == p) && !(p == -1)
            println("WARNING: Using inherent characteristic of the boundary matrix!")
        end
    end

    # Compute the connection matrix in reordered form

    cmRedP = deepcopy(bndA)
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

    # Return the connection matrix information
    
    cm = ConleyMorseCM{typeof(cmMatr.zero)}(
                cmMatr, cmCols, cmPoset, cmLabels, cmMorseSets, cmPoincare)

    if returnbasis
        return cm, basisdict
    else
        return cm
    end
end

