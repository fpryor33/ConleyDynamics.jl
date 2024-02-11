export example_BKMW20_fig1

"""
    lc, mvf = example_BKMW20_fig1()

Create the simplicial complex and multivector field
for the example from Figure 1 in the FoCM 2020 paper
by *Batko, Kaczynski, Mrozek, and Wanner*.

# Examples
```jldoctest
julia> lc, mvf = example_BKMW20_fig1();

julia> cm = connection_matrix(lc, mvf, p=2);

julia> sparse_show(cm.cm)
[0   0   0   0   1]
[0   0   0   0   0]
[0   0   0   0   1]
[0   0   0   0   0]
[0   0   0   0   0]

julia> print(cm.labels)
["A", "AD", "F", "BF", "DE"]
```
"""
function example_BKMW20_fig1()
    # Create the simplicial complex and multivector field for the example
    # from Figure 1 in the FoCM 2020 paper by Batko, Kaczynski, Mrozek,
    # and Wanner.
    
    # Create the simplicial complex

    labels = ["A", "B", "C", "D", "E", "F"]
    intsimplices = [[1,3],[1,4],[3,4],[2,5],[2,6],[5,6],[4,5]]
    lc = create_simplicial_complex(labels, intsimplices)

    # Create the multivector field

    mvf = [["A","AD"], ["D","CD"], ["C","AC"], ["B","BE"], ["E","EF"]]

    # Return the example data

    return lc, mvf
end

