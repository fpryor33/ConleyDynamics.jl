export restrict_dynamics

"""
    restrict_dynamics(lc::LefschetzComplex, mvf::CellSubsets, lcsub::Cells)

Restrict a multivector field to a Lefschetz subcomplex.

For a given multivector field `mvf` on a Lefschetz complex `lc`, and a subcomplex
which is given by the locally closed set represented by `lcsub`, create the
associated Lefschetz subcomplex `lcreduced` and the induced multivector
field `mvfreduced` on the subcomplex. The multivectors of the new multivector
field are the intersections of the original multivectors and the subcomplex.
"""
function restrict_dynamics(lc::LefschetzComplex, mvf::CellSubsets, lcsub::Cells)
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

    # Make sure the cell subset is locally closed

    if !lefschetz_is_locally_closed(lc, lcsubI)
        error("The cell subset has to be locally closed!")
    end

    # Create the restricted multivector field

    mvfredI = Vector{Vector{Int}}()

    for cellsub in mvfI
        cellsubred = intersect(cellsub, lcsubI)
        if length(cellsubred)>1
            push!(mvfredI,cellsubred)
        end
    end

    mvfredL = convert_cellsubsets(lc, mvfredI)

    # Create the reduced Lefschetz complex

    lcred = lefschetz_subcomplex(lc, lcsubI)
 
    # Return the results

    if return_as_labels
        mvfred = mvfredL
    else
        mvfred = convert_cellsubsets(lcred, mvfredL)
    end

    return lcred, mvfred
end

