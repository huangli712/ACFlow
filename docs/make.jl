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
        "Internals" => Any[
            "README" => "internals/README.md",
            "ACFlow APIs" => Any[
                "ACFlow" => "internals/apis/acflow.md",
                "Global" => "internals/apis/global.md",
            ],
        ],
        "Theory" => Any[],
    ],
)
