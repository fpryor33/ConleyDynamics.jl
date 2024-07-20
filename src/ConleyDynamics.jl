"""
    module ConleyDynamics

Collection of tools for computational Conley theory.
"""
module ConleyDynamics

using Graphs
using Combinatorics
using LinearAlgebra
using Random
using Printf
using Luxor
using Colors

# Load the type definitions

include("./sparse/sparse_type.jl")
include("./conley/composite_types.jl")

# Load the methods

include("./conley/cm_reduce.jl")
include("./conley/connection_matrix.jl")
include("./conley/conley_index.jl")
include("./conley/morse_sets.jl")
include("./conley/invariant_hull.jl")
include("./conley/restrict_dynamics.jl")
include("./conley/remove_exit_set.jl")

include("./homology/ph_reduce.jl")
include("./homology/persistent_homology.jl")
include("./homology/homology.jl")
include("./homology/relative_homology.jl")

include("./utils/convert_cells.jl")
include("./utils/convert_coordinates.jl")
include("./utils/create_simplicial_complex.jl")
include("./utils/create_simplicial_rectangle.jl")
include("./utils/create_simplicial_delaunay.jl")
include("./utils/create_cubical_complex.jl")
include("./utils/create_cubical_rectangle.jl")
include("./utils/create_cubical_box.jl")
include("./utils/permute_lefschetz_complex.jl")
include("./utils/lefschetz_boundary.jl")
include("./utils/lefschetz_openhull.jl")
include("./utils/lefschetz_closure.jl")
include("./utils/lefschetz_lchull.jl")
include("./utils/lefschetz_is_closed.jl")
include("./utils/lefschetz_is_locally_closed.jl")
include("./utils/lefschetz_clomo_pair.jl")
include("./utils/lefschetz_skeleton.jl")
include("./utils/lefschetz_subcomplex.jl")
include("./utils/lefschetz_closed_subcomplex.jl")
include("./utils/lefschetz_filtration.jl")
include("./utils/lefschetz_gfp_conversion.jl")
include("./utils/lefschetz_field.jl")
include("./utils/create_mvf_hull.jl")
include("./utils/create_planar_mvf.jl")
include("./utils/create_spatial_mvf.jl")
include("./utils/manifold_boundary.jl")
include("./utils/extract_multivectors.jl")
include("./utils/planar_nontransverse_edges.jl")

include("./sparse/sparse_basic_functions.jl")
include("./sparse/sparse_from_lists.jl")
include("./sparse/sparse_from_full.jl")
include("./sparse/sparse_access.jl")
include("./sparse/sparse_remove.jl")
include("./sparse/sparse_permute.jl")
include("./sparse/sparse_minor.jl")
include("./sparse/sparse_elementary_op.jl")
include("./sparse/sparse_is_sut.jl")
include("./sparse/sparse_multiply.jl")

include("./examples/example_MW_fig01.jl")
include("./examples/example_MW_fig02.jl")
include("./examples/example_MW_fig03.jl")
include("./examples/example_MW_fig04.jl")
include("./examples/example_MW_fig11.jl")
include("./examples/example_BKMW20_fig1.jl")
include("./examples/example_BKMW20_fig3.jl")
include("./examples/example_nonunique.jl")
include("./examples/example_critical_simplex.jl")
include("./examples/example_moebius.jl")

include("./plots/plot_planar_simplicial.jl")
include("./plots/plot_planar_simplicial_morse.jl")
include("./plots/plot_planar_cubical.jl")
include("./plots/plot_planar_cubical_morse.jl")

end

