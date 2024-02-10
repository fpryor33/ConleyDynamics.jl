export lefschetz_skeleton

"""
    lefschetz_skeleton(lc::LefschetzComplex, subcomp::Vector{Int}, skdim::Int)

Compute the `skdim`-dimensional skeleton of a Lefschetz complex subset.

The computed skeleton is for the closure of the subcomplex given by `subcomp`.
"""
function lefschetz_skeleton(lc::LefschetzComplex, subcomp::Vector{Int}, skdim::Int)
    #
    # Compute the skdim-dimensional skeleton  of a Lefschetz complex subset
    #

    # Extract the boundary matrix as a sparse matrix

    if typeof(lc.boundary)==Matrix
        bnd = sparse_from_full(lc.boundary)
    else
        bnd = lc.boundary
    end

    # Compute the closure of the subcomplex
    
    lcclosure = lefschetz_closure(lc, subcomp)

    # Extract the cells of dimension skdim

    lcskeleton = Vector{Int}()

    for k in lcclosure
        if lc.dimensions[k] == skdim
            push!(lcskeleton,k)
        end
    end

    # Return the skeleton indices

    return lcskeleton
end

"""
    lefschetz_skeleton(lc::LefschetzComplex, subcomp::Vector{String}, skdim::Int)

Compute the `skdim`-dimensional skeleton of a Lefschetz complex subset.

The computed skeleton is for the closure of the subcomplex given by `subcomp`.
"""
function lefschetz_skeleton(lc::LefschetzComplex, subcomp::Vector{String}, skdim::Int)
    #
    # Compute the skdim-dimensional skeleton  of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcskeleton = lefschetz_skeleton(lc, subcompI, skdim)
    return lcskeleton
end

