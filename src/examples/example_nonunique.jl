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

julia> cm1 = connection_matrix(lc1, mvf);

julia> cm2 = connection_matrix(lc2, mvf);

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
["B", "G", "GI", "BI", "DE", "FG", "AFH", "CDI", "GHI"]
julia> sparse.show(cm2.cm)
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
["B", "H", "GH", "BI", "DE", "FG", "AFH", "CDI", "GHI"]
```
"""
function example_nonunique()
    # Create two representations of a simplicial complex and one
    # multivector field which illustrates nonunique connection matrices

    # Initialize the simplicial complex in natural order
    
    labels = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    intsimplices = [[1,2,8],[2,8,9],[2,3,9],
                    [1,6,8],[6,7,8],[7,8,9],[4,7,9],[3,4,9],
                    [5,6,7],[4,5,7]]
    strsimplices = convert_simplices(intsimplices, labels)

    # Create the multivector field using string labels

    mvf = [["A","AB"], ["C","BC"], ["D","CD"], ["E","EF"], ["F","AF"],
           ["G","GI"], ["I","HI"], ["H","GH"],
           ["AH","ABH"], ["BH","BHI"], ["CI","BCI"], ["DI","DGI"],
           ["DG","DEG"], ["EG","EFG"], ["FH","FGH"]]

    # Create the first simplicial complex
    
    labels1 = labels
    lc1 = create_simplicial_complex(labels1, strsimplices)

    # Create the second simplicial complex, via vertex reordering

    labels2 = ["A", "B", "C", "D", "E", "F", "H", "I", "G"]
    lc2 = create_simplicial_complex(labels2, strsimplices)

    # Return the example data

    return lc1, lc2, mvf
end

