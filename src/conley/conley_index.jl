export conley_index

# Largest set size for which `algorithm="auto"` uses the rank algorithm.
# See the docstring of `conley_index` for the rationale behind this value.

const CONLEY_INDEX_RANK_MAX = 100

"""
    conley_index(lc::AbstractComplex, subcomp::Vector{Int};
                 [algorithm::String="auto", validate::Bool=true])

Determine the Conley index of a Lefschetz complex subset.

The Conley index of the locally closed set `subcomp` is returned as a
vector `betti` of Betti numbers, where `betti[k]` is the Betti number
in dimension `k-1`. The computations are performed over the field
associated with the Lefschetz complex boundary matrix. Since the
Conley index of a locally closed set is the Lefschetz homology of
the set itself, only the cells in `subcomp` are used, and the
computational cost does not depend on the size of `lc`. The
following optional arguments are supported:

* `algorithm::String` selects the homology algorithm:
  - `algorithm = "auto"` chooses between `rank` and `morse` based
    on the size of `subcomp`, see below. This is the default.
  - `algorithm = "morse"` selects an algorithm based on Morse
    reductions,
  - `algorithm = "pmorse"` selects a parallelized algorithm based
    on Morse reductions,
  - `algorithm = "matrix"` selects the standard matrix algorithm,
    as for example described in the book by Edelsbrunner, Harer,
  - `algorithm = "rank"` determines the Betti numbers from the
    ranks of the boundary maps restricted to `subcomp`, which are
    computed via row reduction. This avoids the construction of a
    subcomplex, and it is the fastest choice for small sets such
    as individual multivectors, while `morse` is much faster
    for large sets.
* `validate::Bool` determines whether the function checks that
  `subcomp` is locally closed, and raises an error otherwise.
  The default is `true`. Passing `validate=false` skips this check,
  which should only be done if it is known that the set is locally
  closed, since otherwise the result is meaningless.

All algorithms return the same Betti numbers, they only differ in
their speed. The choice made by `algorithm = "auto"` is based on the
following observations. The `rank` algorithm row reduces the boundary
matrix restricted to `subcomp`, without any further setup. For small
sets this is very fast, but due to fill-in during the row reduction
its cost grows roughly quadratically with the size of the set. The
`morse` algorithm first extracts the subcomplex and constructs a Morse
matching, which adds a fixed overhead, but afterwards its cost grows
almost linearly. In benchmarks on cubical and simplicial complexes,
`rank` is three to seven times faster than `morse` for sets with up
to about ten cells, which is the typical size of a multivector. The
two algorithms take about the same time for sets with 100 to 150
cells, and for larger sets `morse` is faster by a quickly growing
factor, for example by a factor of 5 to 10 for 1000 cells, and by a
factor of more than 100 for 50000 cells. Therefore `auto` uses `rank`
for sets with at most 100 cells, and `morse` for larger sets. Near
this threshold the two algorithms differ by less than a factor of
1.5, so its precise value is not critical. The decision only uses
the number of cells in `subcomp`, and therefore has no cost.
"""
function conley_index(lc::AbstractComplex, subcomp::Vector{Int};
                      algorithm::String="auto", validate::Bool=true)
    #
    # Determine the Conley index of a Lefschetz complex subset
    #

    # Check the algorithm

    if !(algorithm in ["auto", "morse", "pmorse", "matrix", "rank"])
        error("Unknown algorithm " * algorithm * "!")
    end

    # Choose the algorithm based on the size of the set

    if algorithm == "auto"
        algorithm = length(subcomp) <= CONLEY_INDEX_RANK_MAX ? "rank" : "morse"
    end

    # Deal with the trivial case first

    betti = fill(Int(0), lc.dim+1)
    if length(subcomp) == 0
        return betti
    end

    # Do we have a locally closed set?

    if validate && !lefschetz_is_locally_closed(lc, subcomp)
        error("This is not a locally closed set!")
    end

    # For the standard homology algorithms, compute the
    # Lefschetz homology of the subcomplex given by subcomp

    if !(algorithm == "rank")
        lcsub = lefschetz_subcomplex(lc, subcomp, validate=false)
        bettisub = homology(lcsub, algorithm=algorithm)
        betti[1:length(bettisub)] .= bettisub
        return betti
    end

    # Otherwise, use the ranks of the restricted boundary maps.
    # The restricted boundary matrix only has nonzero blocks
    # mapping k-cells to (k-1)-cells, and the pivot columns of its
    # reduced row echelon form determine the ranks of these blocks.

    cells = sort(subcomp)
    bnd = lefschetz_subboundary(lc, cells)
    sparse_rref!(bnd)

    ncount = fill(Int(0), lc.dim+2)
    rank   = fill(Int(0), lc.dim+2)
    for c in cells
        ncount[lc.dimensions[c]+1] += 1
    end
    for r in 1:bnd.nrow
        if length(bnd.rows[r]) > 0
            pivotcol = bnd.rows[r][1]
            rank[lc.dimensions[cells[pivotcol]]+1] += 1
        end
    end

    # Compute the Betti numbers and return them

    for k = 1:lc.dim+1
        betti[k] = ncount[k] - rank[k] - rank[k+1]
    end

    return betti
end

"""
    conley_index(lc::AbstractComplex, subcomp::Vector{String};
                 [algorithm::String="auto", validate::Bool=true])

Determine the Conley index of a Lefschetz complex subset.

This method expects the set `subcomp` in label form, and accepts
the same optional arguments `algorithm` and `validate` as the
method for subsets in index form.
"""
function conley_index(lc::AbstractComplex, subcomp::Vector{String};
                      algorithm::String="auto", validate::Bool=true)
    #
    # Determine the Conley index of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    return conley_index(lc, subcompI, algorithm=algorithm, validate=validate)
end
