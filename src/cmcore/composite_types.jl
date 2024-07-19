export LefschetzComplex, ConleyMorseCM, CellList, CellListVector

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
* `matrix::SparseMatrix{T}`: Connection matrix
* `columns::Vector{Int}`: Corresponding columns in the boundary matrix
* `poset::Vector{Int}`: Poset indices for the connection matrix columns
* `labels::Vector{String}`: Labels for the connection matrix columns
* `morse::Vector{Vector{String}}`: Vector of Morse sets in original complex
* `conley::Vector{Vector{Int}}`: Vector of Conley indices for the Morse sets
* `complex::LefschetzComplex`: The Conley complex as a Lefschetz complex
"""
struct ConleyMorseCM{T}
    matrix::SparseMatrix{T}
    columns::Vector{Int}
    poset::Vector{Int}
    labels::Vector{String}
    morse::Vector{Vector{String}}
    conley::Vector{Vector{Int}}
    complex::LefschetzComplex
end

"""
    CellList = Union{Vector{Int},Vector{String}}

A list of cells of a Lefschetz complex.

This data type is used to represent subsets of a Lefschetz
complex. It is used for individual isolated invariant sets,
locally closed subsets, and multivectors.
"""
CellList = Union{Vector{Int},Vector{String}}

"""
    CellListVector = Union{Vector{Vector{Int}},Vector{Vector{String}}}

A vector of cell lists.

This data type is used to represent a collection of subsets of
a Lefschetz complex. It is used for Morse decompositions and
for multivector fields.
"""
CellListVector = Union{Vector{Vector{Int}},Vector{Vector{String}}}

