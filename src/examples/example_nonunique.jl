export example_nonunique

"""
    lc1, lc2, mvf = example_nonunique()

Create two representations of a simplicial complex and one multivector
field which illustrates nonunique connection matrices.

The two complexes `lc1` and `lc2` represent the same simplicial
complex, but differ in the ordering of the labels.

# Examples
```jldoctest
julia> lc1, lc2, mvf = example_nonunique();

julia> cm1 = connection_matrix(lc1, mvf, p=2);

julia> cm2 = connection_matrix(lc2, mvf, p=2);

julia> sparse_show(cm1.cm)
[0   0   0   1   0   1   0   0   0]
[0   0   0   1   0   1   0   0   0]
[0   0   0   0   0   0   0   1   1]
[0   0   0   0   0   0   1   1   0]
[0   0   0   0   0   0   0   1   0]
[0   0   0   0   0   0   1   1   0]
[0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0]

julia> print(cm1.labels)
["1", "6", "68", "18", "34", "56", "057", "238", "678"]
julia> sparse_show(cm2.cm)
[0   0   0   1   0   1   0   0   0]
[0   0   0   1   0   1   0   0   0]
[0   0   0   0   0   0   1   0   1]
[0   0   0   0   0   0   1   1   0]
[0   0   0   0   0   0   0   1   0]
[0   0   0   0   0   0   1   1   0]
[0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0]

julia> print(cm2.labels)
["1", "7", "67", "18", "34", "56", "057", "238", "678"]
```
"""
function example_nonunique()
    # Create two representations of a simplicial complex and one
    # multivector field which illustrates nonunique connection matrices

    # Initialize the simplicial complex in natural order
    
    labels = ["0", "1", "2", "3", "4", "5", "6", "7", "8"]
    intsimplices = [[1,2,8],[2,8,9],[2,3,9],
                    [1,6,8],[6,7,8],[7,8,9],[4,7,9],[3,4,9],
                    [5,6,7],[4,5,7]]
    strsimplices = convert_simplices(intsimplices, labels)

    # Create the multivector field using string labels

    mvf = [["0","01"], ["2","12"], ["3","23"], ["4","45"], ["5","05"],
           ["6","68"], ["8","78"], ["7","67"],
           ["07","017"], ["17","178"], ["28","128"], ["38","368"],
           ["36","346"], ["46","456"], ["57","567"]]

    # Create the first simplicial complex
    
    labels1 = labels
    lc1 = create_simplicial_complex(labels1, strsimplices)

    # Create the second simplicial complex, via vertex reordering

    labels2 = ["0", "1", "2", "3", "4", "5", "7", "8", "6"]
    lc2 = create_simplicial_complex(labels2, strsimplices)

    # Return the example data

    return lc1, lc2, mvf
end

