export convert_mvf

"""
    convert_mvf(lc::LefschetzComplex, mvf::Vector{Vector{Int}})

Convert multivector field `mvf` on the Lefschetz complex `lc` from 
index form to label form.
"""
function convert_mvf(lc::LefschetzComplex, mvf::Vector{Vector{Int}})
    #
    # Convert a multivector field from index to label form
    #

    newmvf = Vector{Vector{String}}()

    for k=1:length(mvf)
        push!(newmvf,Vector{String}())
        for m=1:length(mvf[k])
            push!(newmvf[k],lc.labels[mvf[k][m]])
        end
    end

    return newmvf
end

"""
    convert_mvf(lc::LefschetzComplex, mvf::Vector{Vector{String}})

Convert multivector field `mvf` on the Lefschetz complex `lc` from 
label form to index form.
"""
function convert_mvf(lc::LefschetzComplex, mvf::Vector{Vector{String}})
    #
    # Convert a multivector field from label to index form
    #

    newmvf = Vector{Vector{Int}}()

    for k=1:length(mvf)
        push!(newmvf,Vector{Int}())
        for m=1:length(mvf[k])
            push!(newmvf[k],lc.indices[mvf[k][m]])
        end
    end

    return newmvf
end

