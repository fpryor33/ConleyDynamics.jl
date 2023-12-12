export LefschetzComplex, ConleyMorseCM

"""
    LefschetzComplex{Tpoincare}

Collect the Lefschetz complex information in a struct.

The struct has the following fields:
* `ncells::Int``: number of cells
* `boundary::Matrix{Int}`: boundary matrix, columns give the cell boundaries
* `label::Vector{String}`: vector of labels associated with cell indices
* `index::Dict{String,Int}`: dictionary for finding cell index from label
* `poincare::Vector{Tpoincare}`: vector of Poincare polynomials for the cells
"""
struct LefschetzComplex{Tpoincare}
    ncells::Int
    boundary::Matrix{Int}
    label::Vector{String}
    index::Dict{String,Int}
    poincare::Vector{Tpoincare}
end

"""
    ConleyMorseCM{Tmatrix,Tpoincare}

Collect the connection matrix information in a struct.

The struct has the following fields:
* `cm::Tmatrix`: connection matrix
* `columns::Vector{Int}`: corresponding columns in the boundary matrix
* `poset::Vector{Int}`: poset indices for the connection matrix columns
* `labels::Vector{String}`: labels for the connection matrix columns
* `morsesets::Vector{Vector{String}}`: vector of Morse sets in original complex
* `poincare::Vector{Tpoincare}`: vector of Poincare polynomials for the Morse sets
"""
struct ConleyMorseCM{Tmatrix,Tpoincare}
    cm::Tmatrix
    columns::Vector{Int}
    poset::Vector{Int}
    labels::Vector{String}
    morsesets::Vector{Vector{String}}
    poincare::Vector{Tpoincare}
end

