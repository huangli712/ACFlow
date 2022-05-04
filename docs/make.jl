push!(LOAD_PATH, ENV["ZEN_CORE"])

using Documenter, ZenCore

makedocs(
    sitename = "Zen",
    clean = false,
    authors = "Li Huang",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
    ),
    modules = [ZenCore],
    pages = [
        "Home" => "index.md",
        "Introduction" => "intro.md",
        "Getting started" => Any[
            "README" => "start/README.md",
            "Download Zen" => "start/download.md",
            "Compile Zen" => "start/compile.md",
            "Configure Zen" => "start/configure.md",
        ],
        "Tutorials" => Any[
            "README" => "tutor/README.md",
            "SrVO3" => Any[],
            "FeSe" => Any[],
            "NiO" => Any[],
            "Ce" => Any[],
        ],
        "Guide" => Any[
            "README" => "guide/README.md",
            "Core applications" => Any[
                "How to use" => Any[
                    "Interactive mode" => "guide/core/repl.md",
                    "Batch Mode" => "guide/core/batch.md",
                    "Pipeline" => "guide/core/pipeline.md",
                ],    
                "Input parameters" => Any[
                    "PCASE block" => "guide/core/case.md",
                    "PDFT block" => "guide/core/dft.md",
                    "PDMFT block" => "guide/core/dmft.md",
                    "PIMP block" => "guide/core/impurity.md",
                    "PSOLVER block" => "guide/core/solver.md",
                ],    
            ],
            "Auxiliary tools" => Any[],
            "Components" => Any[
                "Density functional theory" => Any[],
                "Dynamical mean-field theory" => Any[],
                "Quantum impurity solver" => Any[],
                "Kohn-Sham Adaptor" => Any[],
                "Intermediate representation" => Any[],
                "Self-energy functions" => Any[],
                "Mixer" => Any[],
            ],
            "Files" => Any[
                "Standard output" => Any[],
                "case.cycle" => Any[],
                "case.log" => Any[],
                "case.stop" => Any[],
                "case.test" => Any[],    
            ],
        ],
        "Internals" => Any[
            "README" => "internals/README.md",
            "Zen's framework" => Any[],
            "ZenCore APIs" => Any[
                "ZenCore" => "internals/apis/zencore.md",
                "Global" => "internals/apis/global.md",
                "Util" => "internals/apis/util.md",
                "Tetra" => "internals/apis/tetra.md",
                "Types" => "internals/apis/types.md",
                "Config" => "internals/apis/config.md",
                "Base" => "internals/apis/base.md",
                "VASP" => "internals/apis/vasp.md",
                "QE" => "internals/apis/qe.md",
                "PLO" => "internals/apis/plo.md",
                "Wannier" => "internals/apis/wannier.md",
                "IR" => "internals/apis/ir.md",
                "DMFT" => "internals/apis/dmft.md",
                "Solver" => "internals/apis/solver.md",
                "Sigma" => "internals/apis/sigma.md",
                "Mixer" => "internals/apis/mixer.md",
            ],
        ],
        "Theory" => Any[],
    ],
)
