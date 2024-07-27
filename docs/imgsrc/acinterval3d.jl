
using ConleyDynamics

function allencahn3d(x::Vector{Float64})
    #
    # Allen-Cahn projection
    #
    lambda = 3.0 * pi
    c      = lambda / pi
    x1, x2, x3 = x
    y1 = (lambda-1)*x1 - 1.5*c * (x1*x1*x1-x1*x1*x3+x2*x2*x3+2*x1*(x2*x2+x3*x3))
    y2 = (lambda-4)*x2 - 1.5*c * x2 * (2*x1*x1+x2*x2+2*x1*x3+2*x3*x3)
    y3 = (lambda-9)*x3 + 0.5*c * (x1*(x1*x1-3*x2*x2)-3*x3*(2*x1*x1+2*x2*x2+x3*x3))
    return [y1, y2, y3]
end

N = 51
bmax = [1.8, 1.5, 1.0]
lc, coordsI = create_cubical_box(N,N,N);
coordsN = convert_spatial_coordinates(coordsI, -bmax, bmax);
mvf = create_spatial_mvf(lc, coordsN, allencahn3d);

morsedecomp = morse_sets(lc, mvf);
cindices = [conley_index(lc, morseset) for morseset in morsedecomp]
iindices = [1,2,3]

morseinterval = morse_interval(lc, mvf, morsedecomp[iindices]);
lci, mvfi = restrict_dynamics(lc, mvf, morseinterval);

cmi = connection_matrix(lci, mvfi);
cmi.conley

# Create graphics

using GLMakie, Colors

csubsets = vcat([convert_cells(lc,morseinterval)],cmi.morse)

vcolors = zeros(UInt8, N, N, N);
for k = 1:length(csubsets)
    for cube in csubsets[k]
        cubeinfo = cube_information(cube)
        if cubeinfo[7] == 3
            vcolors[cubeinfo[1]+1,cubeinfo[2]+1,cubeinfo[3]+1] = k
        end
    end
end

col1 = colorant"royalblue4"
col2 = colorant"royalblue3"
col3 = colorant"steelblue1"
cols = distinguishable_colors(length(csubsets)-1,
                              [col1,col2,col3], dropseed=true)
alph = fill(Float64(1.0),length(csubsets)-1)
cols = vcat(col3,cols)
alph = vcat(0.25, alph)
colalph = [(cols[k],alph[k]) for k in eachindex(cols)]

f, a, p = voxels(-bmax[1]..bmax[1], -bmax[2]..bmax[2],
                 -bmax[3]..bmax[3], vcolors, color = colalph)

# save("/Users/wanner/Desktop/acinterval3d_3_51a.png",f)
# save("/Users/wanner/Desktop/acinterval3d_3_5ab.png",f)

