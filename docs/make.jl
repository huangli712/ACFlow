push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Documenter, ACFlow, Random

makedocs(
    sitename = "ACFlow",
    clean = false,
    authors = "Li Huang",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
    ),
    modules = [ACFlow],
    pages = [
        "Home" => "index.md",
        "Introduction" => "intro.md",
        "Features" => Any[],
        "Usages" => Any[],
        "Tutorials" => Any[],
        "Theory" => Any[],
        "Internals" => Any[
            "README" => "internals/README.md",
            "ACFlow APIs" => Any[
                "ACFlow" => "internals/apis/acflow.md",
                "Constants" => "internals/apis/global.md",
                "Types" => "internals/apis/types.md",
                "Core" => "internals/apis/base.md",
                "Solvers" => "internals/apis/solver.md",
                "Grid" => "internals/apis/grid.md",
                "Mesh" => "internals/apis/mesh.md",
                "Models" => "internals/apis/model.md",
                "Kernels" => "internals/apis/kernel.md",
                "Configuration" => "internals/apis/config.md",
                "Input and output" => "internals/apis/inout.md",
                "Math" => "internals/apis/math.md",
                "Utilities" => "internals/apis/util.md",
            ],
        ],
    ],
)
