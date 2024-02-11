export plot_planar_simplicial

"""
    plot_planar_simplicial(sc::LefschetzComplex,
                           coords::Vector{<:Vector{<:Real}},
                           fname::String;
                           mvf::MultiVectorField,
                           labeldir::Vector{<:Real},
                           hfac::Real,
                           vfac::Real)

Create an svg image of a planar simplicial complex, and if
specified, a Forman vector field on it.

The vector `coords` contains coordinates for every one of the
vertices of the simplicial complex `sc`. The image will be saved
in the file with name `fname`. The optional vector `labeldir`
contains directions for the vertex labels. These have to be
reals between 0 and 4, with 0,1,2,3 corresponding to E,N,W,S.
If `mvf` is specified and is a Forman vector field, then this
Forman vector field is drawn as well. The optional constants
`hfac` and `vfac` contain the horizontal and vertical scale
vectors.
"""
function plot_planar_simplicial(sc::LefschetzComplex,
                                coords::Vector{<:Vector{<:Real}},
                                fname::String;
                                mvf::MultiVectorField=Vector{Vector{Int}}([]),
                                labeldir::Vector{<:Real}=Vector{Int}([]),
                                hfac::Real=1.5,
                                vfac::Real=1.5)
    #
    # Create an svg image of a planar simplicial complex
    #

    # Create proper coordinates
    
    coordsx = sort(unique([coords[k][1] for k in 1:length(coords)]))
    coordsy = sort(unique([-coords[k][2] for k in 1:length(coords)]))
    coordsw = coordsx[end] - coordsx[1]
    coordsh = coordsy[end] - coordsy[1]
    centerx = div(coordsx[end]+coordsx[1], 2)
    centery = div(coordsy[end]+coordsy[1], 2)
    
    pcoords = [[coords[k][1]-centerx, -coords[k][2]-centery]
               for k in 1:length(coords)]
    
    figw = Int(round(coordsw * hfac))
    figh = Int(round(coordsh * vfac))

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
    
    @svg begin
    
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
                      offset=8)
            end
        end

    end figw figh fname
end

