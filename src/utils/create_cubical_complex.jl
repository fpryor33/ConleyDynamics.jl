export create_cubical_complex, create_cube_label

"""
    create_cubical_complex(cubes::Vector{String})

Initialize a Lefschetz complex from a cubical complex.

The vector `cubes` contains a list of all the highest-dimensional cubes
necessary to define the cubical complex. Every cube is represented as a
string as follows:

* d integers, which correspond to the coordinates of a point in d-dimensional Euclidean space
* a point `.`
* d integers 0 or 1, which give the interval length in the respective dimension

The first d integers all have to occupy the same number of characters. In addition,
if the occupied space is L characters for each coordinate, the coordinates only
can take values from 0 to 10^L - 2. This is due to the fact that the boundary
operator will add one to certain coordinates, and they still need to be 
representable withing the same L digits.

For example, the string `030600.101` corresponds to the point `(3,6,0)` in
three dimensions. The dimensions are 1, 0, and 1, and therefore this string
corresponds to the cube `[3,4] x [6] x [0,1]`. The same cube could have also
been represented by `360.101` or by `003006000.101`.

!!! warning
    Note that the labels all have to have the same format!
"""
function create_cubical_complex(cubes::Vector{String})
    #
    # Create a Lefschetz complex struct for a cubical complex.
    #

    # Determine the point dimension and the coordinate field width

    pointdim = length(cubes[1]) - first(findfirst(".",cubes[1]))
    pointlen = Int((length(cubes[1]) - 1) / pointdim - 1)

    # Create a dictionary with all label to integer data information
    # 
    # If label is any label in cubes, then labelintinfodict[label]
    # is an integer vector with the following entries:
    #
    # 1:pointdim:              coordinates of the anchor point
    # 1+pointdim:2*pointdim:   interval length in each dimension
    # 1+2*pointdim:            dimension of the cube

    labelintinfodict = Dict{String,Vector{Int}}()
    for k = 1:length(cubes)
        curlabel = cubes[k]
        curintinfo = Vector{Int}()
        for m = 1:pointdim
            i1 = 1 + (m-1) * pointlen
            i2 = m * pointlen
            push!(curintinfo,parse(Int,curlabel[i1:i2]))
        end
        for m = 1:pointdim
            i1 = pointdim * pointlen + 1 + m
            push!(curintinfo,parse(Int,curlabel[i1]))
        end
        curdim = sum(curintinfo[pointdim+1:2*pointdim])
        push!(curintinfo,curdim)

        if minimum(curintinfo[1:pointdim]) < 0
            error("Point coordinates have to be at least 0!")
        end
        if maximum(curintinfo[1:pointdim]) > 10^pointlen - 2
            error("Point coordinates overflow! (Cannot be 10^w-1)")
        end

        labelintinfodict[curlabel] = curintinfo
    end

    # Find the dimension of the cubical complex
    
    CCdim = maximum([labelintinfodict[n][1+2*pointdim] for n in cubes])

    # Create labels for all cubes in the complex, organized by
    # dimensions, as well as ordered vectors of labels

    labelsets = [Set{String}() for _ in 0:CCdim]
    labelvect = [Vector{String}() for _ in 0:CCdim]
    
    for k = 1:length(cubes)      # Initialize with the given cubes
        curcube = cubes[k]
        curcdim = labelintinfodict[curcube][1+2*pointdim]
        push!(labelsets[curcdim+1], curcube)
    end

    # Work your way from highest to lowest dimension, keep adding faces

    for k = CCdim:-1:1
        
        # Create vector of labels at dimension k

        labelvect[k+1] = sort(collect(labelsets[k+1]))

        # Loop through the cubes and add boundaries

        for curcube in labelvect[k+1]

            curintinfo = labelintinfodict[curcube]
            curones = findall(x -> x==1, curintinfo[1+pointdim:2*pointdim])

            for m=1:k
                newintinfoP = deepcopy(curintinfo)
                newintinfoN = deepcopy(curintinfo)
                ione = curones[m]
                newintinfoP[ione] += 1
                newintinfoP[ione+pointdim] = 0
                newintinfoN[ione+pointdim] = 0
                newintinfoP[1+2*pointdim] -= 1
                newintinfoN[1+2*pointdim] -= 1
                newlabelP = create_cube_label(pointdim, pointlen, newintinfoP)
                newlabelN = create_cube_label(pointdim, pointlen, newintinfoN)
                push!(labelsets[k], newlabelP)
                push!(labelsets[k], newlabelN)
                labelintinfodict[newlabelP] = newintinfoP
                labelintinfodict[newlabelN] = newintinfoN
            end
        end
    end

    labelvect[1] = sort(collect(labelsets[1]))

    # After these preparations, create the Lefschetz complex structure!!

    # Organize cube labels in one vector and create label to index dictionary

    CClabelvec = reduce(vcat,labelvect)
    CClabelindexdict = Dict{String,Int}(CClabelvec[j] => j for j=1:length(CClabelvec))
    CCn = length(CClabelvec)

    # Create the vector of dimensions for the cubes

    CCdimvec = Vector{Int}()
    for k=0:CCdim
        for m=1:length(labelvect[k+1])
            push!(CCdimvec,k)
        end
    end

    # Create the boundary map

    Br = Vector{Int}()
    Bc = Vector{Int}()
    Bv = Vector{Int}()
    for k = 1:CCn
        if CCdimvec[k] > 0

            curcdim = CCdimvec[k]
            curcube = CClabelvec[k]
            curintinfo = labelintinfodict[curcube]
            curones = findall(x -> x==1, curintinfo[1+pointdim:2*pointdim])
            coeff = Int(1)

            for m=1:curcdim
                newintinfoP = deepcopy(curintinfo)
                newintinfoN = deepcopy(curintinfo)
                ione = curones[m]
                newintinfoP[ione] += 1
                newintinfoP[ione+pointdim] = 0
                newintinfoN[ione+pointdim] = 0
                newintinfoP[1+2*pointdim] -= 1
                newintinfoN[1+2*pointdim] -= 1
                newlabelP = create_cube_label(pointdim, pointlen, newintinfoP)
                newlabelN = create_cube_label(pointdim, pointlen, newintinfoN)
                bcellP = CClabelindexdict[newlabelP]
                bcellN = CClabelindexdict[newlabelN]

                push!(Br,bcellP)   # Row index
                push!(Bc,k)        # Column index
                push!(Bv,coeff)    # Matrix entry
                
                push!(Br,bcellN)   # Row index
                push!(Bc,k)        # Column index
                push!(Bv,-coeff)   # Matrix entry
                
                coeff = -coeff
            end
        end
    end

    B = sparse_from_lists(CCn,CCn,0,Int(0),Int(1),Br,Bc,Bv)

    # Create the Lefschetz complex

    lc = LefschetzComplex(CCn,CCdim,B,CClabelvec,CClabelindexdict,CCdimvec)

    # Return the Lefschetz complex

    return lc
end

function create_cube_label(pointdim::Int, pointlen::Int, newintinfo::Vector{Int})
    #
    # Create the label for a cube
    #

    labformat = "%0" * string(pointlen) * "d"
    labelvec = Vector{String}()

    for k = 1:pointdim
        push!(labelvec, Printf.format(Printf.Format(labformat), newintinfo[k]))
    end

    push!(labelvec, ".")

    for k = 1:pointdim
        push!(labelvec, string(newintinfo[k+pointdim]))
    end

    return join(labelvec)
end

