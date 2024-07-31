export convert_cells, convert_cellsubsets

"""
    convert_cells(lc::LefschetzComplex, cl::Vector{Int})

Convert cell list `cl` in the Lefschetz complex `lc` from 
index form to label form.
"""
function convert_cells(lc::LefschetzComplex, cl::Vector{Int})
    #
    # Convert a cell list from index to label form
    #

    return lc.labels[cl]
end

"""
    convert_cells(lc::LefschetzComplex, cl::Vector{String})

Convert cell list `cl` in the Lefschetz complex `lc` from 
label form to index form.
"""
function convert_cells(lc::LefschetzComplex, cl::Vector{String})
    #
    # Convert a cell list from label to index form
    #

    return [lc.indices[k] for k in cl]
end

"""
    convert_cellsubsets(lc::LefschetzComplex, clsub::Vector{Vector{Int}})

Convert CellSubsets `clsub` in the Lefschetz complex `lc` from 
index form to label form.
"""
function convert_cellsubsets(lc::LefschetzComplex, clsub::Vector{Vector{Int}})
    #
    # Convert a CellSubsets from index to label form
    #

    newclsub = Vector{Vector{String}}()

    for k=1:length(clsub)
        push!(newclsub,Vector{String}())
        for m=1:length(clsub[k])
            push!(newclsub[k],lc.labels[clsub[k][m]])
        end
    end

    return newclsub
end

"""
    convert_cellsubsets(lc::LefschetzComplex, clsub::Vector{Vector{String}})

Convert CellSubsets `clsub` in the Lefschetz complex `lc` from 
label form to index form.
"""
function convert_cellsubsets(lc::LefschetzComplex, clsub::Vector{Vector{String}})
    #
    # Convert a CellSubsets from label to index form
    #

    newclsub = Vector{Vector{Int}}()

    for k=1:length(clsub)
        push!(newclsub,Vector{Int}())
        for m=1:length(clsub[k])
            push!(newclsub[k],lc.indices[clsub[k][m]])
        end
    end

    return newclsub
end

