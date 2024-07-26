# The following code implements a 3d projection of the
# Allen-Cahn equation subject to homogeneous Dirichlet
# boundary conditions on the domain (0,pi). The projection
# is onto the span of the first three sine eigenfunctions
# of the Dirichlet Laplacian.
#
# Here are a number of sample runs:
#
# * lambda = 3*pi and N = 25:
#   This takes about 10 seconds and gives the correct
#   Morse decomposition. There are 2 equilibria each of
#   index 0, 1, and 2, and one equilibrium of index 3.
#
# There are also a few other interesting runs:
#
# * lambda = 3*pi and N = 41:
#   This takes 46 seconds and recovers the full Morse
#   decomposition.
#
# * lambda = 3*pi and N = 51:
#   This takes 120 seconds, gives everything.
#
# * lambda = 3*pi and N = 21:
#   This takes 7.5 seconds but cannot resolve the full
#   Morse decomposition. It only gives 6 Morse sets.
#   There are 2 equilibria each of index 0 and 2, and
#   one of index 3. But the two equilibria of index 1
#   are part of one Morse set.
#
# * lambda = 4*pi and N = 51:
#   This run takes about 370 seconds and gives a total
#   of eleven equilibria: 1 index 3, 2 index 2, 4 index 1,
#   and 4 index 0. Four of these are spurious!
#

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

N = 25
bmax = [1.8, 1.5, 1.0]
lc, coordsI = create_cubical_box(N,N,N);
coordsN = convert_spatial_coordinates(coordsI, -bmax, bmax);
mvf = create_spatial_mvf(lc, coordsN, allencahn3d);

morsedecomp = morse_sets(lc, mvf);
morseinterval = morse_interval(lc, mvf, morsedecomp);
lci, mvfi = restrict_dynamics(lc, mvf, morseinterval);

cmi = connection_matrix(lci, mvfi);
cmi.conley

# Create graphics

using GLMakie

colors = zeros(Int64, N, N, N);
for k = 1:length(cmi.morse)
    for cube in cmi.morse[k]
        cubeinfo = cube_information(cube)
        if cubeinfo[7] == 3
            colors[cubeinfo[1]+1,cubeinfo[2]+1,cubeinfo[3]+1] = k
        end
    end
end

f, a, p = voxels(-bmax[1]..bmax[1], -bmax[2]..bmax[2],
                 -bmax[3]..bmax[3], colors,
                 colorrange = (1, length(cmi.morse)),
                 colormap = [:blue, :red],
                 is_air = x -> (x == 0))
save("/Users/wanner/Desktop/allencahn3d_3_25.png",f)

