"""
    module ConnectionMatrices

Computation of connection matrices over finite fields.
"""
module ConnectionMatrices

using Nemo
using Graphs
using Combinatorics
using Random

include("./core/composite_types.jl")
include("./core/convert_matrix_int.jl")
include("./core/convert_matrix_gfp.jl")
include("./core/renumber_poset.jl")
include("./core/admissible_order.jl")
include("./core/update_low.jl")
include("./core/homogeneous_columns.jl")
include("./core/target_columns.jl")
include("./core/cm_columns.jl")
include("./core/cm_create.jl")
include("./core/connection_matrix.jl")
include("./core/create_simplicial_complex.jl")
include("./core/permute_lefschetz_complex.jl")
include("./core/convert_mvf.jl")
include("./core/convert_simplices.jl")

include("./sparse/type_sparse.jl")
include("./sparse/sparse_from_lists.jl")

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

end

