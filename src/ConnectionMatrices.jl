"""
    module ConnectionMatrices

Computation of connection matrices over finite fields.
"""
module ConnectionMatrices

using Nemo
using Graphs
using Combinatorics

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
include("./core/convert_mvf.jl")
include("./core/convert_simplices.jl")

include("./examples/example_MW_fig1.jl")
include("./examples/example_MW_fig2.jl")
include("./examples/example_MW_fig4.jl")
include("./examples/example_BKMW20_fig1.jl")
include("./examples/example_BKMW20_fig3.jl")
include("./examples/example_nonunique.jl")
include("./examples/example_critical_simplex.jl")

end

