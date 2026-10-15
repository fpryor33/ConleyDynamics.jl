# The algorithm behind `morse_sets`

This note explains how `morse_sets` in `src/conley/morse_sets.jl` works. The
first part covers the function as a whole. The second part looks in detail at
how the Hasse diagram of the Morse poset is computed when the function is
called with `poset=true`.

```julia
morse_sets(lc::AbstractComplex, mvf::CellSubsets; poset::Bool=false)
```

The function returns the nontrivial Morse sets of the multivector field `mvf`
on the Lefschetz complex `lc`. With `poset=true` it also returns the adjacency
matrix of the Hasse diagram of the Morse sets.


## Part 1: Overview of the algorithm

### Step 1: Convert to integer form

If `mvf` is given by cell labels, it is converted to cell indices with
`convert_cellsubsets`. The flag `return_as_labels` remembers this, so the
Morse sets can be converted back to labels at the end.

### Step 2: Map cells to multivectors

`cell2mv[c]` is the index of the multivector containing cell `c`. While this
table is built, the function checks that no cell lies in two multivectors. If
one does, it raises the error "The multivector field is not a partition!".

### Step 3: Build the dynamics digraph

The combinatorial flow is encoded in a directed graph `dg` with one vertex per
cell and two kinds of edges:

- **Boundary edges:** for every nonzero entry `(r, c)` of the boundary matrix
  there is an edge `c → r`, from a cell to its face. The edges come straight
  from the sparse storage via `lists_from_sparse`, which costs time linear in
  the number of nonzero entries.
- **Multivector cycles:** the cells of each multivector `[m1, m2, ..., mk]`
  are joined by a directed cycle `m1 → m2 → ... → mk → m1`. This lets the flow
  move freely within a multivector while using only `k` edges instead of `k²`.
  A one-cell multivector gets a self-loop, so every cell in the field is a
  vertex of the graph.

A path in `dg` is exactly a path of the multivector flow: from a cell you can
move to any cell of your own multivector, or to any face of your cell.

### Step 4: Strongly connected components

`strongly_connected_components_tarjan(dg)` splits the cells into strongly
connected components (SCCs). Each SCC is a union of multivectors, because each
multivector lies on a cycle. Each SCC is sorted for a reproducible output. This
step costs time linear in the size of the graph.

### Step 5: Decide which SCCs are nontrivial

An SCC is a Morse set if it carries nontrivial dynamics:

- **Two or more multivectors:** if the SCC contains cells from at least two
  different multivectors (`length(unique(cell2mv[scc[k]])) > 1`), it contains
  a genuine recurrent orbit and is always nontrivial.
- **A single multivector:** it is locally closed, and it is nontrivial exactly
  when its Conley index is nonzero. The index is computed with
  `conley_index(lc, scc[k], algorithm="auto", validate=false)`. This works
  only on the set itself, so the cost depends on the size of the
  multivector, not of the whole complex. The `"auto"` algorithm uses the
  rank-based method for multivectors with at most 100 cells, where it is
  fastest, and the Morse-reduction method for larger ones, where the rank
  method would grow quadratically (see the `conley_index` docstring).
  Validation is skipped because local closedness is guaranteed.

The SCCs are independent of each other, so this loop runs with
`Threads.@threads`. Each iteration writes only its own entry of
`isnontrivial`.

### Step 6: Extract the Morse sets

`ntindices` holds the indices of the nontrivial SCCs, and
`morsesets = scc[ntindices]`. The order of the Morse sets is the order in
which Tarjan's algorithm returns the components.

### Step 7 (optional): Hasse diagram of the Morse poset

See Part 2.

### Step 8: Return

If the input was given by labels, the Morse sets are converted back to labels.
The function returns `morsesets`, or `(morsesets, hasse)` when `poset=true`.

### Overall cost

Without the poset, all steps except the Conley index calls are linear in the
number of cells plus the number of nonzero boundary entries. The index
computations are local to each single-multivector SCC, so in total
`morse_sets` is close to linear in the size of the complex.


## Part 2: Computing the Hasse diagram

### The problem

Morse set M1 lies above Morse set M2 if `dg` has a path from M1 to M2. The
Hasse diagram keeps only the covering pairs: M1 → M2 where no third Morse set
lies strictly between them.

The output is a dense `nms × nms` Boolean matrix `hasse`, where `nms` is the
number of Morse sets. `hasse[m2, m1] = true` means there is a Hasse edge
`m1 → m2`, i.e. M1 covers M2.

### Why the earlier implementation was slow

An earlier version called `has_path(dg, u, v)` for every pair of Morse sets.
Each call is a breadth-first search over the whole cell graph, so this cost
O(nms² · (N + E)) for N cells and E edges. When there are many Morse sets this
is cubic in N. It then ran `transitivereduction` on the dense, transitively
closed graph, which can add another O(nms³). The current implementation
replaces the pairwise searches with a single pass over the graph, done in an
order where each piece of information is ready by the time it is needed. It
produces exactly the same output.

### Step 1: Lookup tables

```julia
cell2scc[v] = k        # cell v lies in SCC k
scc2ms[k]   = m        # SCC k is Morse set m (0 if it is not a Morse set)
```

`ntindices` already lists which SCCs became Morse sets, in the order of
`morsesets`. So `scc2ms[ntindices] = 1:nms` numbers them consistently with the
output.

### Step 2: Processing order

```julia
sccorder = reverse(topological_sort(condensation(dg, scc)))
```

`condensation` collapses each SCC to a single vertex. The result has no
cycles: there is an edge k → k' if some cell in SCC k has an edge into SCC k'.
A topological sort lists every SCC before all its successors, so the reversed
list puts every SCC after all its successors. Processing in this order means
anything needed from a successor has already been computed.

Tarjan's algorithm also returns the components in this order, but the code
does not rely on that implementation detail. Computing the order explicitly
costs only O(N + E).

### Step 3: First-hit sets

For each SCC `k`, `firsthit[k]` is the set of Morse sets you can reach from
`k` along a path that does not pass through any other Morse set on the way.
These are the Morse sets "directly below" `k`, with only trivial SCCs in
between.

```julia
for v in scc[k], w in outneighbors(dg, v)
    kw = cell2scc[w]
    kw == k && continue                    # edge stays inside the SCC
    targets = scc2ms[kw] > 0 ? (scc2ms[kw]:scc2ms[kw]) : firsthit[kw]
    ...add targets to firsthit[k]...
end
```

For each edge leaving the SCC into a successor `kw`, there are two cases:

- **`kw` is a Morse set:** the path stops there, so that Morse set is a first
  hit. We don't look past it, because whatever lies beyond it is recorded in
  its own `firsthit`.
- **`kw` is a trivial SCC:** the path passes through it, so we inherit
  everything it hits first: `firsthit[kw]`. This is already complete because
  of the processing order.

The range `m:m` is just a one-element collection, so both cases share one
loop.

**The `stamp` trick:** one Morse set can be reached through many edges, but
each should appear only once in the list. `stamp[m] == k` means "already added
while processing SCC `k`". Because `k` differs for each SCC, the array never
has to be reset, so removing duplicates costs O(1) per candidate. Neither
sorting nor `unique` is needed.

### Step 4: From first hits to the Hasse diagram

**Claim 1:** following first-hit edges repeatedly gives exactly the order
between Morse sets. Take any path from M1 to M2 and list the Morse sets it
passes through: M1, Q1, Q2, …, M2. Each consecutive pair is a first-hit pair,
because the path between them avoids other Morse sets. Conversely, each
first-hit pair corresponds to a real path.

**Claim 2:** every cover of M1 is in `firsthit(M1)`. Suppose M1 covers M2. On
any path from M1 to M2, the first Morse set reached must be M2 itself.
Otherwise it would be some Q with M1 > Q > M2, contradicting the cover.

So only the first-hit candidates need to be checked. A candidate is discarded
exactly when it can also be reached through another candidate:

```julia
for k in sccorder                      # successors first
    m1 = scc2ms[k]; m1 == 0 && continue
    fill!(covered, false)
    for m2 in firsthit[k]
        covered .|= reach[m2]          # everything strictly below some first hit
    end
    for m2 in firsthit[k]
        reach[m1][m2] = true
        if !covered[m2]                # not reachable via another first hit
            hasse[m2,m1] = true        # => M1 covers M2
        end
    end
    reach[m1] .|= covered              # reach(M1) = firsthit ∪ (their reach)
end
```

- `reach[m]` is a bit vector of all Morse sets strictly below `m`. Because
  Morse sets are processed successors-first, `reach[m2]` is complete by the
  time it is used.
- `covered` collects everything strictly below some first-hit set. If a
  candidate `m2` is in there, then M1 → Q → … → M2 for some other first-hit
  set Q, so the edge M1 → M2 is redundant. A Morse set cannot reach itself
  (the order is antisymmetric), so `m2` can never mark itself as covered.
- If `m2` is not covered, nothing lies between M1 and M2, so it is a Hasse
  edge.

**Small example:** say Morse set A has paths to both B and C that avoid other
Morse sets, and B has a path to C. Then `firsthit(A) = {B, C}`. When A is
processed, `reach[B]` contains C, so C is covered. Only A → B goes into the
Hasse diagram, and C hangs below B through B's own entry.

### Why the output matches the earlier implementation

The Hasse diagram of a partial order is unique, so any correct method gives
the same matrix. The Morse sets and their order are computed before this
block and are not touched by it. The orientation convention
`hasse[m2, m1] = true` for an edge m1 → m2 is the same as before. The test set
`"Morse set poset"` in `test/test_conley.jl` checks the result against an
independent reference, which finds reachability with a plain search and takes
the covers as R minus R∘R.

### Cost

- **Steps 1–2:** O(N + E).
- **Step 3:** O(N + E) times the typical size of a first-hit set, which is
  usually small.
- **Step 4:** one bit-vector OR per first-hit pair, roughly (number of
  first-hit pairs) × nms/64 machine operations.
- **Memory:** nms² bits for `reach`, which is less than the dense `hasse`
  output already takes.

In the extreme test case where every cell of a cubical grid is its own Morse
set (so nms = N), the poset computation for 2601 cells dropped from about
1.7 s to 0.003 s, and 40401 cells take about 0.2 s. What still grows faster
than linearly there comes from the `nms × nms` tables, which the dense return
type requires anyway.
