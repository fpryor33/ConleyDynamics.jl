export example_BKMW20_fig3

"""
    lcomplex, mvf = example_BKMW20_fig3()

Create the simplicial complex and multivector field for the example
from Figure 3 in the FoCM 2020 paper by *Batko, Kaczynski, Mrozek,
and Wanner*.
"""
function example_BKMW20_fig3()
    # Create the simplicial complex and multivector field for the example
    # from Figure 3 in the FoCM 2020 paper by Batko, Kaczynski, Mrozek,
    # and Wanner.
    
    # Create the simplicial complex

    labels = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    intsimplices = [[1,4,5],[1,2,5],[2,5,6],[2,3,6],[3,6,7],
                    [4,5,8],[5,8,9],[5,6,9],[6,9,10],[6,7,10]]
    strsimplices = convert_simplices(intsimplices, labels)
    lc = create_simplicial_complex(labels, strsimplices)

    # Create the multivector field

    mvf = [["A","AD"], ["B","AB"], ["C","BC"], ["F","FG"], ["G","GJ"],
           ["H","DH"], ["I","FI"], ["J","IJ"],
           ["AE","ABE"], ["BE","BEF"], ["CF","BCF"], ["CG","CFG"],
           ["DE","DEH"], ["EH","EHI"], ["EI","EFI"], ["FJ","FIJ"]]

    # Return the example data

    return lc, mvf
end

