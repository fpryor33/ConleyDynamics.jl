export convert_clist, convert_clistvec

"""
    convert_clist(lc::LefschetzComplex, cl::Vector{Int})

Convert cell list `cl` in the Lefschetz complex `lc` from 
index form to label form.
"""
function convert_clist(lc::LefschetzComplex, cl::Vector{Int})
    #
    # Convert a cell list from index to label form
    #

    return lc.labels[cl]
end

"""
    convert_clist(lc::LefschetzComplex, cl::Vector{String})

Convert cell list `cl` in the Lefschetz complex `lc` from 
label form to index form.
"""
function convert_clist(lc::LefschetzComplex, cl::Vector{String})
    #
    # Convert a cell list from label to index form
    #

    return [lc.indices[k] for k in cl]
end

"""
    convert_clistvec(lc::LefschetzComplex, clv::Vector{Vector{Int}})

Convert CellListVector `clv` in the Lefschetz complex `lc` from 
index form to label form.
"""
function convert_clistvec(lc::LefschetzComplex, clv::Vector{Vector{Int}})
    #
    # Convert a CellListVector from index to label form
    #

    newclv = Vector{Vector{String}}()

    for k=1:length(clv)
        push!(newclv,Vector{String}())
        for m=1:length(clv[k])
            push!(newclv[k],lc.labels[clv[k][m]])
        end
    end

    return newclv
end

"""
    convert_clistvec(lc::LefschetzComplex, clv::Vector{Vector{String}})

Convert CellListVector `clv` in the Lefschetz complex `lc` from 
label form to index form.
"""
function convert_clistvec(lc::LefschetzComplex, clv::Vector{Vector{String}})
    #
    # Convert a CellListVector from label to index form
    #

    newclv = Vector{Vector{Int}}()

    for k=1:length(clv)
        push!(newclv,Vector{Int}())
        for m=1:length(clv[k])
            push!(newclv[k],lc.indices[clv[k][m]])
        end
    end

    return newclv
end

