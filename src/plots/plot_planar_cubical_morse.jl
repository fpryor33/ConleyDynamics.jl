export plot_planar_cubical_morse

"""
    plot_planar_cubical_morse(cc::LefschetzComplex,
                              fname::String,
                              morsesets::MultiVectorField;
                              [hfac::Real=1.2,]
                              [vfac::Real=1.2,]
                              [cubefac::Real=0,]
                              [pdim::Vector{Bool}=[false,true,true],]
                              [pv::Bool=false])

Create an image of a planar cubical complex, together with
Morse sets, or also selected multivectors.

The image will be saved in the file with name `fname`, and the ending
determines the image type. Accepted are `.pdf`, `.svg`, `.png`, and `.eps`.

The vector `morsesets` contains a list of Morse sets, or more general,
subsets of the cubical complex. For every `k`, the set described
by `morsesets[k]` will be shown in a distinct color.

The optional constants `hfac` and `vfac` contain the horizontal and
vertical scale vectors for the margins, while `cubefac` describes a uniform
scale. If `cubefac=0` the latter is automatically determined. The vector `pdim`
specifies in which dimensions cells are drawn; the default only shows edges
and squares. Finally if one passes the argument `pv=true`, then in addition
to saving the file a preview is displayed.
"""
function plot_planar_cubical_morse(cc::LefschetzComplex,
                                   fname::String,
                                   morsesets::MultiVectorField;
                                   hfac::Real=1.2,
                                   vfac::Real=1.2,
                                   cubefac::Real=0,
                                   pdim::Vector{Bool}=[false,true,true],
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
        cubefac = 800.0 / maximum([1, oxmax-oxmin, oymax-oymin])
    end

    # Compute the coordinates

    coords = Vector{Vector{Float64}}()
    for k = 1:length(vertices)
        intinfo = cube_information(cc.labels[k])
        push!(coords, [cubefac*intinfo[1], cubefac*intinfo[2]])
    end

    # Create proper coordinates
   
    cx0 = minimum([c[1] for c in coords])
    cx1 = maximum([c[1] for c in coords])
    cy0 = minimum([c[2] for c in coords])
    cy1 = maximum([c[2] for c in coords])
    figw = Int(round((cx1 - cx0) * hfac))
    figh = Int(round((cy1 - cy0) * vfac))
    figdx = (cx1 - cx0) * (hfac-1.0) * 0.5
    figdy = (cy1 - cy0) * (vfac-1.0) * 0.5

    pcoords = [[figdx + (c[1] - cx0) * (figw-2.0*figdx) / (cx1-cx0),
                figdy + (cy1 - c[2]) * (figh-2.0*figdy) / (cy1-cy0)]
               for c in coords]

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
    
    # Plot the cubical complex

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

    # Plot the Morse sets

    if morsesets isa Vector{Vector{Int}}
        msI = morsesets
    else
        msI = convert_mvf(morsesets, cc)
    end

    col1 = colorant"royalblue4"
    col2 = colorant"royalblue3"
    col3 = colorant"steelblue1"
    cols = distinguishable_colors(length(msI), [col1,col2,col3], dropseed=true)

    for m in eachindex(msI)
        setcolor(cols[m])
        for k in msI[m]
            cdim = cc.dimensions[k]
            if (cdim == 0) & pdim[1]
                circle(points[k], 5, action = :fill)
            elseif (cdim == 1) & pdim[2]
                k1 = cellvertices[k][1]
                k2 = cellvertices[k][2]
                line(points[k1],points[k2])
                strokepath()
            elseif (cdim == 2) & pdim[3]
                k1 = cellvertices[k][1]
                k2 = cellvertices[k][3]
                k3 = cellvertices[k][4]
                k4 = cellvertices[k][2]
                poly([points[k1],points[k2],points[k3],points[k4]],
                           action = :fill; close=true)
            end
        end
    end

    # Finish the drawing, and preview if desired

    finish()
    if pv
        preview()
    end
end

