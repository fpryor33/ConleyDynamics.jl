export plot_planar_simplicial_morse

"""
    plot_planar_simplicial_morse(sc::LefschetzComplex,
                                 coords::Vector{<:Vector{<:Real}},
                                 fname::String,
                                 morsesets::CellListVector;
                                 [hfac::Real=1.2,]
                                 [vfac::Real=1.2,]
                                 [sfac::Real=0,]
                                 [pdim::Vector{Bool}=[false,true,true],]
                                 [pv::Bool=false])

Create an image of a planar simplicial complex, together with
Morse sets, or also selected multivectors.

The vector `coords` contains coordinates for every one of the
vertices of the simplicial complex `sc`. The image will be saved
in the file with name `fname`, and the ending determines the image
type. Accepted are `.pdf`, `.svg`, `.png`, and `.eps`.

The vector `morsesets` contains a list of Morse sets, or more general,
subsets of the simplicial complex. For every `k`, the set described
by `morsesets[k]` will be shown in a distinct color.

The optional constants `hfac` and `vfac` contain the horizontal and
vertical scale vectors for the margins, while `sfac` describes a uniform
scale. If `sfac=0` the latter is automatically determined. The vector `pdim`
specifies in which dimensions cells are drawn; the default only shows edges
and triangles. Finally if one passes the argument `pv=true`, then in addition
to saving the file a preview is displayed.
"""
function plot_planar_simplicial_morse(sc::LefschetzComplex,
                                      coords::Vector{<:Vector{<:Real}},
                                      fname::String,
                                      morsesets::CellListVector;
                                      hfac::Real=1.2,
                                      vfac::Real=1.2,
                                      sfac::Real=0,
                                      pdim::Vector{Bool}=[false,true,true],
                                      pv::Bool=false)
    #
    # Create an image of a planar simplicial complex
    #

    # Create proper coordinates
   
    cx0 = minimum([c[1] for c in coords])
    cx1 = maximum([c[1] for c in coords])
    cy0 = minimum([c[2] for c in coords])
    cy1 = maximum([c[2] for c in coords])

    if iszero(sfac)
        sfac = 800.0 / maximum([1, cx1-cx0, cy1-cy0])
    end

    figw  = Int(round((cx1 - cx0) * hfac * sfac))
    figh  = Int(round((cy1 - cy0) * vfac * sfac))
    figdx = (cx1 - cx0) * (hfac-1.0) * 0.5 * sfac
    figdy = (cy1 - cy0) * (vfac-1.0) * 0.5 * sfac

    pcoords = [[figdx + (c[1] - cx0) * (figw-2.0*figdx) / (cx1-cx0),
                figdy + (cy1 - c[2]) * (figh-2.0*figdy) / (cy1-cy0)]
               for c in coords]

    # Create a list of vertex points

    points = [Point(pcoords[k][1],pcoords[k][2]) for k in 1:length(pcoords)]

    # Create a list of vertices for each simplex

    cellvertices = Vector{Vector{Int}}()
    for k = 1:sc.ncells
        cv = lefschetz_skeleton(sc, [k], 0)
        if !(length(cv) == 1+sc.dimensions[k])
            error("The complex is not a simplicial complex!")
        end
        push!(cellvertices, cv)
    end

    # Create the image
   
    Drawing(figw, figh, fname)
    background("white")
    sethue("black")
    
    # Plot the simplicial complex

    for k = sc.ncells:-1:1
        cdim = sc.dimensions[k]
        if (cdim == 0) && pdim[1]
            setcolor("royalblue4")
            circle(points[k], 5, action = :fill)
        elseif (cdim == 1) && pdim[2]
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            setcolor("royalblue3")
            line(points[k1],points[k2])
            strokepath()
        elseif (cdim == 2) && pdim[3]
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            k3 = cellvertices[k][3]
            setcolor("steelblue1")
            poly([points[k1],points[k2],points[k3]],
                       action = :fill; close=true)
        end
    end

    # Plot the Morse sets

    if morsesets isa Vector{Vector{Int}}
        msI = morsesets
    else
        msI = convert_clistvec(sc, morsesets)
    end

    col1 = colorant"royalblue4"
    col2 = colorant"royalblue3"
    col3 = colorant"steelblue1"
    cols = distinguishable_colors(length(msI), [col1,col2,col3], dropseed=true)

    for m in eachindex(msI)
        setcolor(cols[m])
        setopacity(0.6)
        for k in msI[m]
            cdim = sc.dimensions[k]
            if (cdim == 0) && pdim[1]
                circle(points[k], 5, action = :fill)
            elseif (cdim == 1) && pdim[2]
                k1 = cellvertices[k][1]
                k2 = cellvertices[k][2]
                line(points[k1],points[k2])
                strokepath()
            elseif (cdim == 2) && pdim[3]
                k1 = cellvertices[k][1]
                k2 = cellvertices[k][2]
                k3 = cellvertices[k][3]
                poly([points[k1],points[k2],points[k3]],
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

