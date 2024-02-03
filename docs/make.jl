#
# To create the html documentation locally, use
# the following command:
# 
#    julia --color=yes --project make.jl --local-html
#
# For the Latex pdf documentation, use instead:
#
#    julia --color=yes --project make.jl --latex-pdf
#

push!(LOAD_PATH,"../src/")

using Documenter
using ConleyDynamics

pageslist = ["Overview" => "index.md",
             "Manual" => Any[
                 "man/guide.md",
                 "man/lefschetz.md",
                 "man/connections.md",
                 "man/homology.md",
                 "man/sparse.md",
                 "man/examples.md"
                  ],
             "Core API" => Any[
                 "apicore/datastruct.md",
                 "apicore/utils.md",
                 "apicore/cmcore.md",
                 "apicore/homology.md",
                 "apicore/sparse.md"
                  ],
            ]

DocMeta.setdocmeta!(ConleyDynamics, :DocTestSetup,
                    :(using ConleyDynamics); recursive=true)

if "--local-html" in ARGS
     makedocs(sitename="ConleyDynamics.jl",
        modules=[ConleyDynamics],
        format = Documenter.HTML(prettyurls = false),
        pages = pageslist,
        authors = "Thomas Wanner"
        )
elseif "--latex-pdf" in ARGS
    makedocs(sitename="ConleyDynamics.jl",
        modules=[ConleyDynamics],
        format = Documenter.LaTeX(platform = "docker"),
        pages = pageslist,
        authors = "Thomas Wanner"
        )
else
    makedocs(sitename="ConleyDynamics.jl",
        modules=[ConleyDynamics],
        format = Documenter.HTML(),
        pages = pageslist,
        authors = "Thomas Wanner"
        )
    deploydocs(
        repo = "github.com/almost6heads/ConleyDynamics.jl.git",
        )
end

