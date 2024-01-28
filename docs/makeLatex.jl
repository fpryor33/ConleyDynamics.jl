#
# Start this with the following command:
#    julia --color=yes --project makeLatex.jl
#

push!(LOAD_PATH,"../src/")

using DocumenterLaTeX
using Documenter
using ConleyDynamics

makedocs(format = LaTeX(platform = "docker"),
         authors = "Thomas Wanner")

