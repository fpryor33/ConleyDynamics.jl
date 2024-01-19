export convert_mvf

"""
    newmvf = convert_mvf(mvf::Vector{Vector{Int}}, lc::LefschetzComplex)

Convert multivector field `mvf` on the Lefschetz complex `lc` from 
index form to label form.
"""
function convert_mvf(mvf::Vector{Vector{Int}}, lc::LefschetzComplex)
    #
    # Convert a multivector field from index to label form
    #

    newmvf = Vector{Vector{String}}()

    for k=1:length(mvf)
        push!(newmvf,Vector{String}())
        for m=1:length(mvf[k])
            push!(newmvf[k],lc.label[mvf[k][m]])
        end
    end

    return newmvf
end

"""
    newmvf = convert_mvf(mvf::Vector{Vector{String}}, lc::LefschetzComplex)

Convert multivector field `mvf` on the Lefschetz complex `lc` from 
label form to index form.
"""
function convert_mvf(mvf::Vector{Vector{String}}, lc::LefschetzComplex)
    #
    # Convert a multivector field from label to index form
    #

    newmvf = Vector{Vector{Int}}()

    for k=1:length(mvf)
        push!(newmvf,Vector{Int}())
        for m=1:length(mvf[k])
            push!(newmvf[k],lc.index[mvf[k][m]])
        end
    end

    return newmvf
end

