#
# Start this with the following command:
#    julia --color=yes --project make.jl
#

push!(LOAD_PATH,"../src/")

# using Documenter, ConleyDynamics
# For pdf file generation use instead:
using DocumenterLaTeX, Documenter, ConleyDynamics

makedocs(sitename="ConleyDynamics.jl",
         # format = Documenter.HTML(prettyurls = false),
         # For pdf file generation use instead:
         format = LaTeX(platform = "docker"),
         pages = [
                  "Home" => "index.md",
                  "Manual" => Any[
                      "Guide" => "guide.md",
                      "Examples" => "examples.md",
                      ],
                 ],
         authors = "Thomas Wanner"
        )

