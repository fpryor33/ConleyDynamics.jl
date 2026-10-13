export restrict_dynamics

"""
    restrict_dynamics(lc::AbstractComplex, mvf::CellSubsets, lcsub::Cells;
                      [validate::Bool=true])

Restrict a multivector field to a Lefschetz subcomplex.

For a given multivector field `mvf` on a Lefschetz complex `lc`, and a subcomplex
which is given by the locally closed set represented by `lcsub`, create the
associated Lefschetz subcomplex `lcreduced` and the induced multivector
field `mvfreduced` on the subcomplex. The multivectors of the new multivector
field are the intersections of the original multivectors and the subcomplex.

By default, the function raises an error if `lcsub` is not locally closed.
This check can be skipped by passing the optional argument `validate=false`,
which should only be done if it is already known that the set is locally
closed, for example if it was computed using [`morse_interval`](@ref).
"""
function restrict_dynamics(lc::AbstractComplex, mvf::CellSubsets, lcsub::Cells;
                           validate::Bool=true)
    #
    # Restrict a multivector field to a Lefschetz subcomplex
    #

    # Convert the multivector field and cell subset to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = mvf
        return_as_labels = false
    else
        mvfI = convert_cellsubsets(lc, mvf)
        return_as_labels = true
    end

    if lcsub isa Vector{Int}
        lcsubI = lcsub
    else
        lcsubI = convert_cells(lc, lcsub)
    end

    # Create the reduced Lefschetz complex

    lcred = lefschetz_subcomplex(lc, lcsubI, validate=validate)

    # Map the cells of the subcomplex to their indices in the
    # reduced Lefschetz complex, and all other cells to zero

    newindex = fill(Int(0), lc.ncells)
    for (k, cell) in enumerate(sort(lcsubI))
        newindex[cell] = k
    end

    # Create the restricted multivector field in the
    # index form of the reduced Lefschetz complex

    mvfred = Vector{Vector{Int}}()
    for mv in mvfI
        mvred = Vector{Int}()
        for cell in mv
            if newindex[cell] > 0
                push!(mvred, newindex[cell])
            end
        end
        if length(mvred) > 1
            push!(mvfred, mvred)
        end
    end
 
    # Return the results

    if return_as_labels
        return lcred, convert_cellsubsets(lcred, mvfred)
    else
        return lcred, mvfred
    end
end

