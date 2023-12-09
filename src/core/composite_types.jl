export LefschetzComplex

"""
    LefschetzComplex{Tlabel,Tpoincare}

Collect the Lefschetz complex information in a convenient struct.

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

