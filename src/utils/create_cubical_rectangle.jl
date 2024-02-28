export create_cubical_rectangle

"""
    create_cubical_rectangle(nx::Int, ny::Int)

Create a cubical complex covering a rectangle in the plane.

The rectangle is given by the subset `[0,nx] x [0,ny]` of the plane,
and each unit square gives a two-dimensional cube in the resulting
cubical complex.
"""
function create_cubical_rectangle(nx::Int, ny::Int)
    #
    # Create a Lefschetz complex struct for a cubical rectangle.
    #

    # Make sure that we have at least width one in each direction

    if (nx < 1) || (ny < 1)
        error("Width and height have to be at least 1!")
    end

    # Create the vector of two-dimensional cubes

    pointdim = 2
    pointlen = Int(ceil(log(maximum([nx,ny]) + 2) / log(10)))
    cubes = Vector{String}()

    for k = 0:nx-1
        for m = 0:ny-1
            pointinfo = [k, m, 1, 1]
            clabel = cube_label(pointdim, pointlen, pointinfo)
            push!(cubes, clabel)
        end
    end

    # Create and return the Lefschetz complex

    lc = create_cubical_complex(cubes)

    return lc
end

