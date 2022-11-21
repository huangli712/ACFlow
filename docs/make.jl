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
            "Running Modes" => "man/run.md",
            "Input Files" => "man/input.md",
            "Output Files" => "man/output.md",
            "Parameters" => "man/param.md",
        ],
        "Examples" => Any[
            #"Matsubara Self-Energy Function" => "examples/sigma.md",
            "Matsubara Green's Function" => "examples/green1.md",
            "Imaginary-Time Green's Function" => "examples/green2.md",
            "Current-Current Correlation Function" => "examples/current.md",
        ],
        "Theory" => Any[
            "Basic Principles" => "theory/basic.md",
            "Maximum Entropy Method" => "theory/maxent.md",
            "Stochastic Analytical Continuation 1" => "theory/sac1.md",
            "Stochastic Analytical Continuation 2" => "theory/sac2.md",
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
