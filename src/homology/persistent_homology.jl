export persistent_homology

"""
    phsingles, phpairs = persistent_homology(lc::LefschetzComplex,
                                             filtration::Vector{Int};
                                             p::Int=-1)

Complete the persistent homology of a Lefschetz complex filtration.

The function assumes that the order given by the filtration values
is admissible, i.e., the permuted boundary matrix is strictly
upper triangular. The persistence computation is performed over the
finite field `GF(p)` (for prime `p`) or over the rationals (for `p=0`).
The function returns the starting filtration values for infinite length
persistence intervals in `phsingles`, and the birth- and death-filtration
values for finite length persistence intervals in `phpairs`. If the
Lefschetz complex boundary matrix already has been specialized to a field,
then the optional argument `p` can be omitted.
"""
function persistent_homology(lc::LefschetzComplex, filtration::Vector{Int};
                             p::Int=-1)
    #
    # Compute the persistent homology of a Lefschetz complex filtration
    #

    # Find an admissible permutation of the cells

    fvals = sort(unique(filtration))
    adperm = Vector{Int}([])
    for kv in fvals
        append!(adperm, findall(t -> t==kv, filtration))
    end

    # Create the permuted boundary matrix and make sure it
    # is strictly upper triangular

    bndperm = sparse_permute(lc.boundary, adperm, adperm)

    if !sparse_is_sut(bndperm)
        error("Filtration error!")
        return
    end

    # Convert to finite field and compute the persistence

    if (bndperm isa SparseMatrix{Int}) && (bndperm.char == 0)
        if (p == -1)
            error("Homology over Z is not supported, specify p!")
        elseif (p >= 0)
            bndpermgfp = convert_matrix_gfp(bndperm, p)
        else
            error("Wrong characteristic p!")
        end
    else
        bndpermgfp = bndperm
        if !(bndperm.char == p) && !(p == -1)
            println("WARNING: Using inherent characteristic of the complex boundary!")
        end
    end

    permsingles, permpairs = ph_reduce!(bndpermgfp)

    # Extract the correct persistence intervals based on
    # the original order

    phsingles = [Vector{Int}() for _ in 0:lc.dim]
    phpairs = [Vector{Tuple{Int,Int}}() for _ in 0:lc.dim]

    for k=1:length(permsingles)
        singleindex  = adperm[permsingles[k]]
        singlefilter = filtration[singleindex]
        singledim    = lc.dimensions[singleindex]
        push!(phsingles[1+singledim],singlefilter)
    end

    for k=1:length(permpairs)
        pairindex1  = adperm[permpairs[k][1]]
        pairindex2  = adperm[permpairs[k][2]]
        pairfilter1 = filtration[pairindex1]
        pairfilter2 = filtration[pairindex2]
        if !(pairfilter1 == pairfilter2)
            pairdim = lc.dimensions[pairindex1]
            push!(phpairs[1+pairdim],(pairfilter1,pairfilter2))
        end
    end

    # Return the results

    return phsingles, phpairs
end

