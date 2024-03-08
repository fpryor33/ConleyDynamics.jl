export plot_planar_simplicial

"""
    plot_planar_simplicial(sc::LefschetzComplex,
                           coords::Vector{<:Vector{<:Real}},
                           fname::String;
                           [mvf::MultiVectorField,]
                           [labeldir::Vector{<:Real},]
                           [labeldis::Real,]
                           [hfac::Real,]
                           [vfac::Real,]
                           [pv::Bool])

Create an image of a planar simplicial complex, and if
specified, a Forman vector field on it.

The vector `coords` contains coordinates for every one of the
vertices of the simplicial complex `sc`. The image will be saved
in the file with name `fname`, and the ending determines the image
type. Accepted are `.pdf`, `.svg`, `.png`, and `.eps`.

If the optional `mvf` is specified and is a Forman vector field,
then this Forman vector field is drawn as well. The optional
vector `labeldir` contains directions for the vertex labels,
and `labeldis` the distance from the vertex. The directions
have to be reals between 0 and 4, with 0,1,2,3 corresponding
to E,N,W,S. The optional constants `hfac` and `vfac` contain
the horizontal and vertical scale vectors. Finally if one passes
the argument `pv=true`, then in addition to saving the file
a preview is displayed.

# Examples

Suppose we have created a simplicial complex using the commands

```julia
sc, coords = create_simplicial_delaunay(300, 300, 30, 20)
fname = "sc_plot_test.pdf"
```

Then the following code creates an image of the simplicial complex
without labels, but with a preview:

```julia
plot_planar_simplicial(sc, coords, fname, pv=true)
```

If we want to see the labels, we can use

```julia
ldir = fill(0.5, sc.ncells);
plot_planar_simplicial(sc, coords, fname, labeldir=ldir, labeldis=10, pv=true)
```

This command puts all labels in the North-East direction at a distance of 10.
"""
function plot_planar_simplicial(sc::LefschetzComplex,
                                coords::Vector{<:Vector{<:Real}},
                                fname::String;
                                mvf::MultiVectorField=Vector{Vector{Int}}([]),
                                labeldir::Vector{<:Real}=Vector{Int}([]),
                                labeldis::Real=8,
                                hfac::Real=1.5,
                                vfac::Real=1.5,
                                pv::Bool=false)
    #
    # Create an image of a planar simplicial complex
    #

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

    # Create a list of vertices for each simplex

    cellvertices = Vector{Vector{Int}}()
    for k = 1:sc.ncells
        cv = lefschetz_skeleton(sc, [k], 0)
        if !(length(cv) == 1+sc.dimensions[k])
            error("The complex is not a simplicial complex!")
        end
        push!(cellvertices, cv)
    end

    # If necessary prepare for the Forman vector field plot

    if length(mvf) > 0

        # First make sure that the multivector field has type Int

        if typeof(mvf[1][1]) == String
            mvfI = convert_mvf(mvf, sc)
        else
            mvfI = mvf
        end

        # Find the barycenters

        barycs = Vector{Point}()
        for k = 1:sc.ncells
            cdim = sc.dimensions[k]
            if cdim == 0
                push!(barycs, points[k])
            elseif cdim == 1
                k1 = cellvertices[k][1]
                k2 = cellvertices[k][2]
                push!(barycs, midpoint(points[k1],points[k2]))
            elseif cdim == 2
                k1 = cellvertices[k][1]
                k2 = cellvertices[k][2]
                k3 = cellvertices[k][3]
                push!(barycs, trianglecenter(points[k1],points[k2],points[k3]))
            end
        end

        # Preprocess the Forman vector field

        invector = fill(false,sc.ncells)
        for k = 1:length(mvfI)
            if !(length(mvfI[k]) == 2)
                error("The multivector field is not a Forman vector field!")
            end
            invector[mvfI[k][1]] = true
            invector[mvfI[k][2]] = true
        end
        
        critical = Vector{Int}()
        for k = 1:sc.ncells
            if invector[k] == false
                push!(critical,k)
            end
        end
    end
    
    # Create the image
   
    Drawing(figw, figh, fname)
    background("white")
    sethue("black")
    
    # Plot the simplicial complex

    for k = sc.ncells:-1:1
        cdim = sc.dimensions[k]
        if cdim == 0
            setcolor("royalblue4")
            circle(points[k], 5, action = :fill)
        elseif cdim == 1
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            setcolor("royalblue3")
            line(points[k1],points[k2])
            strokepath()
        elseif cdim == 2
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            k3 = cellvertices[k][3]
            setcolor("steelblue1")
            poly([points[k1],points[k2],points[k3]],
                       action = :fill; close=true)
        end
    end

    if length(mvf) > 0

        # Plot the critical cells

        for m = 1:length(critical)
            k = critical[m]
            setcolor("red3")
            circle(barycs[k], 5, action = :fill)
        end

        # Plot the arrows

        for m = 1:length(mvf)
            k1 = mvfI[m][1]
            k2 = mvfI[m][2]
            cdim1 = sc.dimensions[k1]
            cdim2 = sc.dimensions[k2]
            setcolor("red3")
            if cdim1 < cdim2
                arrow(barycs[k1],barycs[k2],linewidth=3,arrowheadlength=10)
            elseif cdim1 > cdim2
                arrow(barycs[k2],barycs[k1],linewidth=3,arrowheadlength=10)
            else
                error("The multivector field is not a Forman vector field!")
            end
        end
    end

    # Display the vertex labels

    if length(labeldir) > 0
        setcolor("black")
        for k = 1:length(coords)
            label(sc.labels[k], pi*(4.0-labeldir[k])/2.0, points[k],
                  offset=labeldis)
        end
    end

    # Finish the drawing, and preview if desired

    finish()
    if pv
        preview()
    end
end

