export convert_planar_coordinates

"""
    convert_planar_coordinates(coords::Vector{Vector{Float64}},
                               p0::Vector{Float64},
                               p1::Vector{Float64})

Convert a given collection of planar coordinates.

The vector `coords` contains pairs of coordinates, which are then
transformed to fit into the box with vertices `p0 = (p0x,p0y)`
and `p1 = (p1x,p1y)`. It is assumed that `p0` denotes the lower
left box corner, while `p1` is the upper right corner. The function
shifts and scales the coordinates in such a way that every side of
the box contains at least one point. Upon completion, it returns
a new coordinate vector `coordsNew`.

More precisely, if the x-coordinates are spanning the interval
`[xmin,xmax]` and the y-coordinates span `[ymin,ymax]`, then the
point `(x,y)` is transformed to `(xn,yn)` with:

* `xn = p0x + (p1x-p0x) * (x-cxmin) / (cxmax-cxmin)`
* `yn = p0y + (p1y-p0y) * (y-cymin) / (cymax-cymin)`
"""
function convert_planar_coordinates(coords::Vector{Vector{Float64}},
                                    p0::Vector{Float64},
                                    p1::Vector{Float64})
    #
    # Convert planar coordinates to fit inside a given box
    #

    # Find the minimal and maximal coordinates

    cxmin = minimum([c[1] for c in coords])
    cxmax = maximum([c[1] for c in coords])
    cymin = minimum([c[2] for c in coords])
    cymax = maximum([c[2] for c in coords])

    # Create and return the new coordinate vector

    newcoords = Vector{Vector{Float64}}()

    for c in coords
        cnx = p0[1] + (p1[1] - p0[1]) * (c[1] - cxmin) / (cxmax - cxmin)
        cny = p0[2] + (p1[2] - p0[2]) * (c[2] - cymin) / (cymax - cymin)
        push!(newcoords, [cnx,cny])
    end

    return newcoords
end

