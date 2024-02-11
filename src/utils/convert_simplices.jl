export convert_simplices

"""
    convert_simplices(simplices::Vector{Vector{Int}},
                      labels::Vector{String})

Convert list of simplices from index form to label form.
"""
function convert_simplices(simplices::Vector{Vector{Int}},
                           labels::Vector{String})
    #
    # Convert list of simplices from index form to label form
    #
    
    newsimplices = Vector{Vector{String}}()

    for k=1:length(simplices)
        push!(newsimplices,Vector{String}())
        for m=1:length(simplices[k])
            push!(newsimplices[k],labels[simplices[k][m]])
        end
    end

    return newsimplices
end

"""
    convert_simplices(simplices::Vector{Vector{String}},
                      labels::Vector{String})

Convert list of simplices from label form to index form.
"""
function convert_simplices(simplices::Vector{Vector{String}},
                           labels::Vector{String})
    #
    # Convert list of simplices from label form to index form
    #
    
    indices = Dict{String,Int}(labels[k] => k for k in 1:length(labels))
    newsimplices = Vector{Vector{Int}}()

    for k=1:length(simplices)
        push!(newsimplices,Vector{Int}())
        for m=1:length(simplices[k])
            push!(newsimplices[k],indices[simplices[k][m]])
        end
    end

    return newsimplices
end

