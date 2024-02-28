export plot_planar_cubical

"""
    plot_planar_cubical(cc::LefschetzComplex,
                        fname::String;
                        hfac::Real=1.5,
                        vfac::Real=1.5,
                        cubefac::Real=0,
                        pdim::Vector{Bool}=[true,true,true],
                        pv::Bool=false)

Create an image of a planar cubical complex.

The image will be saved in the file with name `fname`, and the
ending determines the image type. Accepted are `.pdf`, `.svg`,
`.png`, and `.eps`. The optional constants `hfac` and `vfac` contain
the horizontal and vertical scale vectors. The optional argument
`cubefac` specifies the side length of an elementary cube for
plotting, and it will be automatically determined otherwise. The
vector `pdim` specifies which cell dimensions should be plotted,
with `pdim[k]` representing dimension `k-1`. Finally if one passes
the argument `pv=true`, then in addition to saving the file
a preview is displayed.
"""
function plot_planar_cubical(cc::LefschetzComplex,
                             fname::String;
                             hfac::Real=1.5,
                             vfac::Real=1.5,
                             cubefac::Real=0,
                             pdim::Vector{Bool}=[true,true,true],
                             pv::Bool=false)
    #
    # Create an image of a planar cubical complex
    #

    # Extract the vertex information

    vertices = lefschetz_skeleton(cc, 0)

    if !(length(vertices) == maximum(vertices))
        error("The vertices are not at the beginning of the cell list!")
    end

    # If desired, determine cubefac automatically

    if iszero(cubefac)
        oxmin = 10^10
        oxmax = 0
        oymin = 10^10
        oymax = 0
        for k = 1:length(vertices)
            intinfo = cube_information(cc.labels[k])
            oxmin = minimum([oxmin, intinfo[1]])
            oxmax = maximum([oxmax, intinfo[1]])
            oymin = minimum([oymin, intinfo[2]])
            oymax = maximum([oymax, intinfo[2]])
        end
        cubefac = 300.0 / maximum([1, oxmax-oxmin, oymax-oymin])
    end

    # Compute the coordinates

    coords = Vector{Vector{Float64}}()
    for k = 1:length(vertices)
        intinfo = cube_information(cc.labels[k])
        push!(coords, [cubefac*intinfo[1], cubefac*intinfo[2]])
    end

    # Create proper coordinates
   
    cxmin = minimum([coords[k][1] for k in 1:length(coords)])
    cxmax = maximum([coords[k][1] for k in 1:length(coords)])
    cymin = minimum([-coords[k][2] for k in 1:length(coords)])
    cymax = maximum([-coords[k][2] for k in 1:length(coords)])
    figw = Int(round((cxmax - cxmin) * hfac))
    figh = Int(round((cymax - cymin) * vfac))
    figwoff = Int(round((cxmax - cxmin) * (hfac-1.0) * 0.5))
    fighoff = Int(round((cymax - cymin) * (vfac-1.0) * 0.5))

    pcoords = [[figwoff + coords[k][1], figh - fighoff - coords[k][2]]
               for k in 1:length(coords)]

    # Create a list of vertex points

    points = [Point(pcoords[k][1],pcoords[k][2]) for k in 1:length(pcoords)]

    # Create a list of vertices for each cube

    cellvertices = Vector{Vector{Int}}()
    for k = 1:cc.ncells
        cv = lefschetz_skeleton(cc, [k], 0)
        if !(length(cv) == 2^cc.dimensions[k])
            error("The complex is not a cubical complex!")
        end
        push!(cellvertices, cv)
    end

    # Create the image
   
    Drawing(figw, figh, fname)
    background("white")
    sethue("black")
    
    # Plot the simplicial complex

    for k = cc.ncells:-1:1
        cdim = cc.dimensions[k]
        if (cdim == 0) & pdim[1]
            setcolor("royalblue4")
            circle(points[k], 5, action = :fill)
        elseif (cdim == 1) & pdim[2]
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            setcolor("royalblue3")
            line(points[k1],points[k2])
            strokepath()
        elseif (cdim == 2) & pdim[3]
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][3]
            k3 = cellvertices[k][4]
            k4 = cellvertices[k][2]
            setcolor("steelblue1")
            poly([points[k1],points[k2],points[k3],points[k4]],
                       action = :fill; close=true)
        end
    end

    # Finish the drawing, and preview if desired

    finish()
    if pv
        preview()
    end
end

