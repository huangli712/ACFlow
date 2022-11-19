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
        "Manual" => Any[
            "Main Features" => "man/feature.md",
            "Implementations" => "man/impl.md",
            "Installation" => "man/install.md",
            "Input Files" => "man/input.md",
            "Output Files" => "man/output.md",
            "Parameters" => "man/param.md",
        ],
        "Examples" => Any[
            "Matsubara Self-Energy Function" => "examples/sigma.md",
            "Imaginary-Time Green's Function" => "examples/green.md",
            "Optical Conductivity" => "examples/optic.md",
        ],
        "Theory" => Any[
            "Basic Principles" => "theory/basic.md",
            "Maximum Entropy Method" => "theory/maxent.md",
            "Stochastic Analytical Continuation" => "theory/sac.md",
            "Stochastic Optimization Method" => "theory/som.md",
        ],
        "Library" => Any[
            "Outline" => "library/outline.md",
            "ACFlow" => "library/acflow.md",
            "Constants" => "library/global.md",
            "Types" => "library/type.md",
            "Core" => "library/base.md",
            "Solvers" => "library/solver.md",
            "Grids" => "library/grid.md",
            "Meshes" => "library/mesh.md",
            "Models" => "library/model.md",
            "Kernels" => "library/kernel.md",
            "Configuration" => "library/config.md",
            "Input and output" => "library/inout.md",
            "Math" => "library/math.md",
            "Utilities" => "library/util.md",
        ],
    ],
)
