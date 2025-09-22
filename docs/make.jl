haskey(ENV,"ACFLOW_HOME") && pushfirst!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Documenter
using Random
using ACFlow

makedocs(
    sitename = "ACFlow: The User Guide",
    clean = true,
    authors = "Li Huang <huangli@caep.cn> and contributors",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
        repolink = "https://github.com/huangli712/ACFlow",
        size_threshold = 409600, # 400kb
        assets = ["assets/acflow.css"],
        collapselevel = 1,
    ),
    #format = Documenter.LaTeX(platform = "none"),
    remotes = nothing,
    modules = [ACFlow],
    pages = [
        "Welcome" => "index.md",
        "Introduction" => Any[
            "Analytic Continuation" => "intro/acp.md",
            "Motivation" => "intro/motivation.md",
            "Acknowledgements" => "intro/ack.md",
            "Citation" => "intro/cite.md",
        ],
        "Manual" => Any[
            "Main Features" => "man/feature.md",
            "Implementations" => "man/impl.md",
            "Installation" => "man/install.md",
            "Running Modes" => "man/run.md",
            "Input Files" => "man/input.md",
            "Output Files" => "man/output.md",
            "Parameters" => "man/param.md",
            "Tricks and Tips" => "man/tricks.md",
            "Graphic User Interface" => "man/gui.md",
            "Benchmark Tools" => "man/test.md",
            "Performance" => "man/performance.md",
        ],
        "Examples" => Any[
            "Matsubara Self-Energy Function" => "examples/sigma.md",
            "Matsubara Green's Function" => "examples/green1.md",
            "Imaginary Time Green's Function" => "examples/green2.md",
            "Current-Current Correlation Function" => "examples/current.md",
        ],
        "Theory" => Any[
            "Basic Principles" => "theory/basic.md",
            "Maximum Entropy Method" => "theory/maxent.md",
            "Barycentric Rational Function" => "theory/rfa.md",
            "Nevanlinna Analytical Continuation" => "theory/nac.md",
            "Stochastic Analytic Continuation 1" => "theory/sac1.md",
            "Stochastic Analytic Continuation 2" => "theory/sac2.md",
            "Stochastic Optimization Method" => "theory/som.md",
            "Stochastic Pole Expansion" => "theory/spx.md",
            "References" => "theory/reference.md",
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
            "Input and Output" => "library/inout.md",
            "Math" => "library/math.md",
            "Utilities" => "library/util.md",
        ],
    ],
)
