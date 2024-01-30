#
# Start this with the following command:
#    julia --color=yes --project make.jl
#

push!(LOAD_PATH,"../src/")

using Documenter
using ConleyDynamics

makedocs(sitename="ConleyDynamics.jl",
         # format = Documenter.HTML(prettyurls = false),
         # For pdf file generation use instead:
         format = Documenter.LaTeX(platform = "docker"),
         pages = [
                  "Home" => "index.md",
                  "Manual" => Any[
                      "Guide" =>               "man/guide.md",
                      "Examples" =>            "man/examples.md",
                      "Data Structures" =>     "man/datastruct.md",
                      "Connection Matrices" => "man/connections.md",
                      "Homology" =>            "man/homology.md",
                      "Sparse Matrices" =>     "man/sparse.md"
                      ],
                 ],
         authors = "Thomas Wanner"
        )

