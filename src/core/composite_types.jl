export LefschetzComplex, ConleyMorseCM

"""
    LefschetzComplex{Tlabel,Tpoincare}

Collect the Lefschetz complex information in a struct.

The struct has the following fields:
* `ncells::Int``: number of cells
* `boundary::Matrix{Int}`: boundary matrix, columns give the cell boundaries
* `label::Vector{Tlabel}`: vector of labels associated with cell indices
* `index::Dict{Tlabel,Int}`: dictionary for finding cell index from label
* `poincare::Vector{Tpoincare}`: vector of Poincare polynomials for the cells
"""
struct LefschetzComplex{Tlabel,Tpoincare}
    ncells::Int
    boundary::Matrix{Int}
    label::Vector{Tlabel}
    index::Dict{Tlabel,Int}
    poincare::Vector{Tpoincare}
end

"""
    ConleyMorseCM{Tmatrix,Tlabel,Tpoincare}

Collect the connection matrix information in a struct.

The struct has the following fields:
* `cm::Tmatrix`: connection matrix
* `columns::Vector{Int}`: corresponding columns in the boundary matrix
* `poset::Vector{Int}`: poset indices for the connection matrix columns
* `labels::Vector{Tlabel}`: labels for the connection matrix columns
* `morsesets::Vector{Vector{Tlabel}}`: vector of Morse sets in original complex
* `poincare::Vector{Tpoincare}`: vector of Poincare polynomials for the Morse sets
"""
struct ConleyMorseCM{Tmatrix,Tlabel,Tpoincare}
    cm::Tmatrix
    columns::Vector{Int}
    poset::Vector{Int}
    labels::Vector{Tlabel}
    morsesets::Vector{Vector{Tlabel}}
    poincare::Vector{Tpoincare}
end

