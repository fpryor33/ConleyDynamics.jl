#
# Start this with the following command:
#    julia --color=yes --project make.jl
#

push!(LOAD_PATH,"../src/")

using Documenter
using ConleyDynamics

makedocs(sitename="ConleyDynamics.jl",
         format = Documenter.HTML(prettyurls = false),
         # For pdf file generation use instead:
         # format = Documenter.LaTeX(platform = "docker"),
         pages = [
                  "Overview" => "index.md",
                  "Manual" => Any[
                      "man/guide.md",
                      "man/examples.md",
                      "man/datastruct.md",
                      "man/connections.md",
                      "man/homology.md",
                      "man/sparse.md"
                      ],
                  "Core API" => Any[
                      "apicore/utils.md",
                      "apicore/cmcore.md",
                      "apicore/homology.md",
                      "apicore/sparse.md"
                     ],
                 ],
         authors = "Thomas Wanner"
        )

