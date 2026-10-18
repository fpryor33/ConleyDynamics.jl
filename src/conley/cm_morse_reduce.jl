export cm_reduce_hms21
export cm_reduce_pmorse26

"""
    cm_reduce_hms21(matrix::SparseMatrix, psetvec::Vector{Int})

Compute the connection matrix.

This function uses the algorithm from the paper Harker, Mischaikow,
Spendlove (Journal of Applied and Computational Topology, 2021).
Assumes that `matrix` is upper triangular and filtered according
to `psetvec`.

# Return values:
* `cmatrix`: Connection matrix
* `cmatrix_cols`: Columns of the connection matrix in the boundary
"""
function cm_reduce_hms21(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Compute the connection matrix
    #

    # Extract the zero and one elements

    tchar = matrix.char
    tzero = matrix.zero
    tone  = matrix.one

    # Find the initial matching

    current_matrix  = matrix
    current_psetvec = psetvec
    current_cvector = Vector{Int}(1:matrix.nrow)

    current_m = matching(current_matrix, current_psetvec)
    
    # If this matching produced fewer critical cells than the
    # number we started out with, we can determine the smaller
    # boundary matrix

    while length(current_m.critical) < current_matrix.nrow

        # Determine the new matrix, whose columns are indexed
        # by the critical cells

        new_matrix, new_cvector = morse_delta(current_matrix, current_m)

        current_matrix  = new_matrix
        current_psetvec = current_psetvec[new_cvector]
        current_cvector = current_cvector[new_cvector]

        # Try to find another matching

        current_m = matching(current_matrix, current_psetvec)
    end

    # Return the result

    return current_matrix, current_cvector
end

"""
    cm_reduce_pmorse26(matrix::SparseMatrix, psetvec::Vector{Int})

Compute the connection matrix.

This function uses a parallelized Morse matching algorithm.
Assumes that `matrix` is upper triangular and filtered according
to `psetvec`.

# Return values:
* `cmatrix`: Connection matrix
* `cmatrix_cols`: Columns of the connection matrix in the boundary
"""
function cm_reduce_pmorse26(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Compute the connection matrix
    #

    # Extract the zero and one elements

    tchar = matrix.char
    tzero = matrix.zero
    tone  = matrix.one

    # Find the initial matching

    current_matrix  = matrix
    current_psetvec = psetvec
    current_cvector = Vector{Int}(1:matrix.nrow)

    current_m = parallel_matching(current_matrix, current_psetvec)
    
    # If this matching produced fewer critical cells than the
    # number we started out with, we can determine the smaller
    # boundary matrix

    while length(current_m.critical) < current_matrix.nrow

        # Determine the new matrix, whose columns are indexed
        # by the critical cells

        new_matrix, new_cvector = parallel_morse_delta(current_matrix, current_m)

        current_matrix  = new_matrix
        current_psetvec = current_psetvec[new_cvector]
        current_cvector = current_cvector[new_cvector]

        # Try to find another matching

        current_m = parallel_matching(current_matrix, current_psetvec)
    end

    # Return the result

    return current_matrix, current_cvector
end

######################################################################################

struct Matching
    critical::Set{Int}
    sources::Set{Int}
    targets::Set{Int}
    arrows::Dict{Int,Int}
end


function matching(matrix::SparseMatrix, filter::Vector{Int})
    #
    # Find a matching for a given boundary matrix, based
    # on a provided filter for the cells.
    #

    @assert filter==sort(filter) "The filter has to be increasing!"
    @assert length(filter)==matrix.nrow "Wrong filter size!"
    @assert matrix.ncol==matrix.nrow "Wrong matrix size!"

    # Initialize the return variables

    critical = Set{Int}()
    sources  = Set{Int}()
    targets  = Set{Int}()
    arrows   = Dict{Int,Int}()

    # Loop through the filter values and create the matching

    for findices in filter_ranges(filter)
        # Find a partial matching for the subcomplex

        cp, ap = partial_matching(matrix, findices)
        
        # Update the return variables

        sp = Vector{Int}()
        tp = Vector{Int}()
        for (s,t) in ap
            push!(sp, s)
            push!(tp, t)
        end

        union!(critical, cp)
        union!(sources, sp)
        union!(targets, tp)
        merge!(arrows, ap)
    end

    return Matching(critical, sources, targets, arrows)
end


function parallel_matching(matrix::SparseMatrix, filter::Vector{Int})
    #
    # Find a matching for a given boundary matrix, based
    # on a provided filter for the cells. This function is
    # parallelized.
    #

    @assert filter==sort(filter) "The filter has to be increasing!"
    @assert length(filter)==matrix.nrow "Wrong filter size!"
    @assert matrix.ncol==matrix.nrow "Wrong matrix size!"

    # Initialize the return variables

    critical = Set{Int}()
    sources  = Set{Int}()
    targets  = Set{Int}()
    arrows   = Dict{Int,Int}()

    # Loop through the filter values and create the matching

    critical_channel = Channel{Vector{Int}}(Inf)
    arrow_channel = Channel{Dict{Int,Int}}(Inf)

    Threads.@threads for findices in filter_ranges(filter)

        # Find a partial matching for the subcomplex

        cp, ap = partial_matching(matrix, findices)
        put!(critical_channel, cp)
        put!(arrow_channel, ap)
    end

    close(critical_channel)
    close(arrow_channel)
    critical_chvec = collect(critical_channel)
    arrow_chvec = collect(arrow_channel)

    # Update the return variables

    for cp in critical_chvec
        union!(critical, cp)
    end

    for ap in arrow_chvec
        sp = Vector{Int}()
        tp = Vector{Int}()
        for (s,t) in ap
            push!(sp, s)
            push!(tp, t)
        end
        union!(sources, sp)
        union!(targets, tp)
        merge!(arrows, ap)
    end

    # Return the results

    return Matching(critical, sources, targets, arrows)
end


function filter_ranges(filter::Vector{Int})
    #
    # Since the filter is increasing, the cells with the same
    # filter value form a range of consecutive indices. Return
    # these index ranges in increasing order of the values.
    #
    franges = Vector{Vector{Int}}()
    n = length(filter)
    i = 1
    while i <= n
        j = i
        while (j < n) && (filter[j+1] == filter[i])
            j += 1
        end
        push!(franges, collect(i:j))
        i = j + 1
    end
    return franges
end


function partial_matching(matrix::SparseMatrix, mind::Vector{Int})
    #
    #  Find a partial matching in the minor `matrix` determined
    #  by the columns in `mind`
    #

    @assert length(mind)>0 "I need something to work with!"
    critical = Vector{Int}([])
    arrows   = Dict{Int,Int}()

    # Create a boundary dictionary for faster manipulation
    # from the given sparse matrix. The dictionary only
    # contains the cells in `mind`, and for every cell it
    # points to its boundary within `mind`.

    mindset = Set{Int}(mind)
    bnddict = Dict{Int, Set{Int}}(c => Set{Int}() for c in mindset)
    for k in mind
        if length(matrix.columns[k]) > 0
            bndset = intersect(Set(matrix.columns[k]), mindset)
            if length(bndset) > 0
                bnddict[k] = bndset
            end
        end
    end
    for k in keys(bnddict)
        if length(bnddict[k]) == 0
            delete!(bnddict, k)
        end
    end

    # Create a second dictionary which records the coboundaries
    # of cells
    
    cbddict = Dict{Int, Set{Int}}(c => Set{Int}() for c in mindset)
    for (k, bndk) in bnddict
        for m in bndk
            push!(cbddict[m], k)
        end
    end
    for k in keys(cbddict)
        if length(cbddict[k]) == 0
            delete!(cbddict, k)
        end
    end

    tomatch = mindset
    removed = Set{Int}()

    # Collect all cells with boundary of size 1

    queue = Vector{Int}()
    for k in keys(bnddict)
        if length(bnddict[k]) == 1
            push!(queue, k)
        end
    end

    # A critical cell is the first cell, in the iteration order of
    # tomatch, without remaining boundary. Since cells are only ever
    # deleted from tomatch, and deleting from a Set does not change
    # the order of the remaining elements, this order is fixed. The
    # cells whose boundary is empty are therefore kept in a heap,
    # keyed by their position in this order. The heap may contain
    # cells which have been removed in the meantime, and these are
    # skipped when they are popped.

    order = collect(tomatch)
    pos   = Dict{Int,Int}(order[i] => i for i in 1:length(order))
    zheap = Vector{Int}()   # max-heap of negated positions
    for k in order
        if !haskey(bnddict, k)
            heap_push!(zheap, -pos[k])
        end
    end

    # Find the matching

    while length(tomatch) > 0

        # Find all coreduction pairs

        while !isempty(queue)
            t = popfirst!(queue)
            
            # Skip if t was already removed, as target of someone else
            (t in removed || length(get(bnddict, t, Set())) != 1) && continue
    
            # Identify the source
            s = first(bnddict[t])
            arrows[s] = t
            
            # Now we have to remove s and t
            for c in (t, s)
                c in removed && error("s error")
                push!(removed, c)
                delete!(tomatch, c)
                
                # Update the boundary of every cell which has c in
                # its boundary
                if haskey(cbddict, c)
                    for cc in cbddict[c]
                        cc in removed && continue

                        # Remove c from the boundary of bc
                        delete!(bnddict[cc], c)
                        
                        # If cc now has a boundary of size 1,
                        # add it to the queue, and if it has
                        # no boundary left, add it to the heap
                        if length(bnddict[cc]) == 1
                            push!(queue, cc)
                        elseif length(bnddict[cc]) == 0
                            heap_push!(zheap, -pos[cc])
                        end
                    end
                end
                
                # Clean up the boundary dictionary entry for
                # the removed cell
                delete!(bnddict, c)
            end
        end

        # Find one element with zero boundary

        if length(tomatch) > 0
            foundit = false
            while !isempty(zheap)
                k = order[-heap_pop!(zheap)]

                # Skip cells which were removed in the meantime
                (k in removed || !isempty(get(bnddict, k, Set{Int}()))) && continue

                # We found a critical cell!
                push!(critical, k)
                foundit = true

                # Update the boundary of every cell which has k in
                # its boundary
                if haskey(cbddict, k)
                    for ck in cbddict[k]
                        ck in removed && continue

                        # Remove k from the boundary of ck
                        delete!(bnddict[ck], k)

                        # If ck now has a boundary of size 1,
                        # add it to the queue, and if it has
                        # no boundary left, add it to the heap
                        if length(bnddict[ck]) == 1
                            push!(queue, ck)
                        elseif length(bnddict[ck]) == 0
                            heap_push!(zheap, -pos[ck])
                        end
                    end
                end

                # Update the cell lists and coboundary dictionary
                push!(removed, k)
                delete!(tomatch, k)
                delete!(bnddict, k)
                break
            end
            @assert foundit "There should have been a critical cell!"
        end
    end

    return critical, arrows
end


struct MorseWorkspace{T}
    x::Vector{T}          # dense work vector
    touched::Vector{Bool} # rows of x that have been modified
    rlist::Vector{Int}    # list of the modified rows
    inheap::Vector{Bool}  # sources which are currently in the heap
    heap::Vector{Int}     # max-heap of sources with nonzero x entry
end

MorseWorkspace(bnd::SparseMatrix) =
    MorseWorkspace(fill(bnd.zero, bnd.nrow), fill(false, bnd.nrow), Vector{Int}(),
                   fill(false, bnd.nrow), Vector{Int}())


function heap_push!(heap::Vector{Int}, v::Int)
    #
    # Insert v into the max-heap
    #
    push!(heap, v)
    i = length(heap)
    while i > 1
        ip = i >> 1
        heap[ip] >= heap[i] && break
        heap[ip], heap[i] = heap[i], heap[ip]
        i = ip
    end
end


function heap_pop!(heap::Vector{Int})
    #
    # Remove and return the largest element of the max-heap
    #
    vmax = heap[1]
    vlast = pop!(heap)
    n = length(heap)
    if n > 0
        heap[1] = vlast
        i = 1
        while true
            il = 2*i
            il > n && break
            ic = (il < n && heap[il+1] > heap[il]) ? il+1 : il
            heap[i] >= heap[ic] && break
            heap[i], heap[ic] = heap[ic], heap[i]
            i = ic
        end
    end
    return vmax
end


function morse_add!(ws::MorseWorkspace, issource::Vector{Bool}, r::Int, val, p::Int)
    #
    # Add val to the work vector entry r, and keep track of the
    # modified rows and the sources with nonzero entry
    #
    ws.x[r] = scalar_add(ws.x[r], val, p)
    if !ws.touched[r]
        ws.touched[r] = true
        push!(ws.rlist, r)
    end
    if issource[r] && !ws.inheap[r] && !iszero(ws.x[r])
        ws.inheap[r] = true
        heap_push!(ws.heap, r)
    end
end


function morse_psi!(ws::MorseWorkspace, bnd::SparseMatrix, issource::Vector{Bool},
                    arrows::Dict{Int,Int}, k::Int)
    #
    # This function implements the Gamma algorithm from the paper
    # by Harker, Mischaikow, Mrozek, Nanda (FOCM, 2014). Note that the
    # version in the paper by Harker, Mischaikow, Spendlove (JACT, 2021)
    # is not correct for general fields.
    #
    # Starting with x = delta e_k, the algorithm repeatedly picks the
    # largest source Q with x[Q] nonzero, and adds the multiple omega
    # of the column delta e_K, where K is matched with Q, which makes
    # x[Q] vanish. If gamma denotes the sum of the multiples omega e_K,
    # then the final vector x equals
    #
    #    delta e_k + delta(gamma(delta e_k)),
    #
    # which is the column needed in morse_delta. On return, the
    # nonzero entries of x are among the rows in ws.rlist. The
    # sources with nonzero entry are kept in a max-heap, which
    # may contain sources whose entry has since become zero.
    #
    p = bnd.char

    for (r, v) in zip(bnd.columns[k], bnd.entries[k])
        morse_add!(ws, issource, r, v, p)
    end

    while !isempty(ws.heap)
        Q = heap_pop!(ws.heap)
        ws.inheap[Q] = false
        iszero(ws.x[Q]) && continue

        K     = arrows[Q]
        omega = scalar_multiply(-ws.x[Q], scalar_inverse(bnd[Q,K],p), p)

        for (jr, vr) in zip(bnd.columns[K], bnd.entries[K])
            morse_add!(ws, issource, jr, scalar_multiply(omega, vr, p), p)
        end
    end
end


function morse_delta_columns!(rlist::Vector{Int}, clist::Vector{Int}, vlist::Vector,
                              ws::MorseWorkspace, bnd::SparseMatrix,
                              issource::Vector{Bool}, arrows::Dict{Int,Int},
                              cproj::Vector{Int}, critpos::Vector{Int}, ells)
    #
    # Compute the columns ells of the reduced boundary matrix, and
    # append their nonzero entries in the critical part to the lists
    #
    for ell in ells
        morse_psi!(ws, bnd, issource, arrows, cproj[ell])
        for r in ws.rlist
            if critpos[r] > 0 && !iszero(ws.x[r])
                push!(rlist, critpos[r])
                push!(clist, ell)
                push!(vlist, ws.x[r])
            end
            ws.x[r]       = bnd.zero
            ws.touched[r] = false
        end
        empty!(ws.rlist)
    end
end


function morse_delta(bnd::SparseMatrix, match::Matching; parallel::Bool=false)
    #
    # This function computes the reduced boundary matrix
    # induced by the specified matching. In essence, we
    # compute for each unit vector v corresponding to a
    # critical cell the vector
    #
    #    psi(delta(phi(x))),
    #
    # as defined in (5) on page 166 in Harker, Mischaikow,
    # Mrozek, Nanda (FOCM, 2014). Notice that their formula
    # can be simplified using the identity
    #
    #    gamma delta gamma = -gamma.
    #
    # Thus, we actually compute
    #
    #    delta v + delta(gamma(delta v)),
    #
    # where delta denotes the boundary operator. This vector
    # is obtained directly in the Gamma algorithm, see morse_psi!.
    # The work per column is proportional to the number of
    # entries visited, not to the size of the matrix.
    #
    @assert bnd.nrow==bnd.ncol "Wrong matrix dimensions!"

    # Initialization

    cproj = sort(collect(match.critical))
    lproj = length(cproj)

    critpos = fill(Int(0), bnd.nrow)
    critpos[cproj] = 1:lproj

    issource = fill(false, bnd.nrow)
    for s in match.sources
        issource[s] = true
    end

    # Compute the columns, either sequentially or split into
    # chunks which each have their own workspace

    nchunks = parallel ? min(Threads.nthreads(), lproj) : min(1, lproj)
    chunks  = [m:nchunks:lproj for m in 1:nchunks]
    rlists  = [Vector{Int}() for _ in 1:nchunks]
    clists  = [Vector{Int}() for _ in 1:nchunks]
    vlists  = [Vector{typeof(bnd.zero)}() for _ in 1:nchunks]

    if parallel
        Threads.@threads for ic in 1:nchunks
            morse_delta_columns!(rlists[ic], clists[ic], vlists[ic], MorseWorkspace(bnd),
                                 bnd, issource, match.arrows, cproj, critpos, chunks[ic])
        end
    else
        for ic in 1:nchunks
            morse_delta_columns!(rlists[ic], clists[ic], vlists[ic], MorseWorkspace(bnd),
                                 bnd, issource, match.arrows, cproj, critpos, chunks[ic])
        end
    end

    # Construct the new boundary matrix

    D = sparse_from_lists(lproj, lproj, bnd.char, bnd.zero, bnd.one,
                          reduce(vcat, rlists; init=Int[]),
                          reduce(vcat, clists; init=Int[]),
                          reduce(vcat, vlists; init=Vector{typeof(bnd.zero)}()))
    return D, cproj
end


parallel_morse_delta(bnd::SparseMatrix, match::Matching) =
    morse_delta(bnd, match; parallel=true)
