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
        "Getting started" => "start.md",
        "Tutorials" => "tutor.md",
        "User's guide" => "guide.md",
        "Theory" => "theory.md",
        "Internals" => Any[
            "README" => "internals/README.md",
            "ACFlow" => "internals/acflow.md",
            "Constants" => "internals/global.md",
            "Types" => "internals/types.md",
            "Core" => "internals/base.md",
            "Solvers" => "internals/solver.md",
            "Grid" => "internals/grid.md",
            "Mesh" => "internals/mesh.md",
            "Models" => "internals/model.md",
            "Kernels" => "internals/kernel.md",
            "Configuration" => "internals/config.md",
            "Input and output" => "internals/inout.md",
            "Math" => "internals/math.md",
            "Utilities" => "internals/util.md",
        ],
    ],
)
