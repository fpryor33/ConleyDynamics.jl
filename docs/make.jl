#
# Start this with the following command:
#    julia --color=yes --project make.jl
#

push!(LOAD_PATH,"../src/")

using Documenter, ConleyDynamics

# makedocs(sitename="ConleyDynamics.jl")
makedocs(sitename="ConleyDynamics.jl", format = Documenter.HTML(prettyurls = false))

