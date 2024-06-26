export LefschetzComplex, ConleyMorseCM, MultiVectorField

"""
    LefschetzComplex

Collect the Lefschetz complex information in a struct.

The struct has the following fields:
* `ncells::Int`: Number of cells
* `dim::Int`: Dimension of the complex
* `boundary::SparseMatrix`: Boundary matrix, columns give the cell boundaries
* `labels::Vector{String}`: Vector of labels associated with cell indices
* `indices::Dict{String,Int}`: Dictionary for finding cell index from label
* `dimensions::Vector{Int}`: Vector cell dimensions
The coefficient field is specified by the boundary matrix.
"""
struct LefschetzComplex
    ncells::Int
    dim::Int
    boundary::SparseMatrix
    labels::Vector{String}
    indices::Dict{String,Int}
    dimensions::Vector{Int}
end

"""
    ConleyMorseCM{T}

Collect the connection matrix information in a struct.

The struct has the following fields:
* `cm::SparseMatrix{T}`: Connection matrix
* `columns::Vector{Int}`: Corresponding columns in the boundary matrix
* `poset::Vector{Int}`: Poset indices for the connection matrix columns
* `labels::Vector{String}`: Labels for the connection matrix columns
* `morsesets::Vector{Vector{String}}`: Vector of Morse sets in original complex
* `poincare::Vector{Vector{Int}}`: Vector of Poincare polynomial coefficients
  for the Morse sets. The k-th entry is the coefficient of t^(k-1).
* `complex::LefschetzComplex`: The Conley complex as a Lefschetz complex.
"""
struct ConleyMorseCM{T}
    cm::SparseMatrix{T}
    columns::Vector{Int}
    poset::Vector{Int}
    labels::Vector{String}
    morsesets::Vector{Vector{String}}
    poincare::Vector{Vector{Int}}
    complex::LefschetzComplex
end

"""
    MultiVectorField = Union{Vector{Vector{Int}},Vector{Vector{String}}}

Type of a multivector field.
"""
MultiVectorField = Union{Vector{Vector{Int}},Vector{Vector{String}}}

