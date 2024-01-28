export LefschetzComplex, ConleyMorseCM, MultiVectorField

"""
    LefschetzComplex

Collect the Lefschetz complex information in a struct.

The struct has the following fields:
* `ncells::Int`: number of cells
* `dim::Int`: dimension of the complex
* `boundary::Matrix{Int}`: boundary matrix, columns give the cell boundaries
* `labels::Vector{String}`: vector of labels associated with cell indices
* `indices::Dict{String,Int}`: dictionary for finding cell index from label
* `dimensions::Vector{Int}`: vector cell dimensions
"""
struct LefschetzComplex
    ncells::Int
    dim::Int
    boundary::Union{Matrix{Int},SparseMatrix{Int}}
    labels::Vector{String}
    indices::Dict{String,Int}
    dimensions::Vector{Int}
end

"""
    ConleyMorseCM{Tmatrix}

Collect the connection matrix information in a struct.

The struct has the following fields:
* `cm::Tmatrix`: connection matrix
* `columns::Vector{Int}`: corresponding columns in the boundary matrix
* `poset::Vector{Int}`: poset indices for the connection matrix columns
* `labels::Vector{String}`: labels for the connection matrix columns
* `morsesets::Vector{Vector{String}}`: vector of Morse sets in original complex
* `poincare::Vector{Vector{Int}}`: vector of Poincare polynomial coefficients
  for the Morse sets. The k-th entry is the coefficient of t^(k-1).
"""
struct ConleyMorseCM{Tmatrix}
    cm::Tmatrix
    columns::Vector{Int}
    poset::Vector{Int}
    labels::Vector{String}
    morsesets::Vector{Vector{String}}
    poincare::Vector{Vector{Int}}
end

"""
    MultiVectorField = Union{Vector{Vector{Int}},Vector{Vector{String}}}

Type of a multivector field.
"""
MultiVectorField = Union{Vector{Vector{Int}},Vector{Vector{String}}}

