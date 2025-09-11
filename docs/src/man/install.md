It is an easy task to install the ACFlow toolkit.

## [Quick Installation](@id install)

First, since it is written in pure Julia language, it is necessary to install the Julia runtime environment in advance. The newest version of Julia is always preferred (version > 1.60). Because the core codes only rely on Julia's built-in standard library, no the third-party packages are needed.

Second, just download source codes of the ACFlow toolkit from its github repository. It should be a compressed file, such as `acflow.zip` or `acflow.tar.gz`. Please uncompress it in your favorite directory by using the following commands:

```shell
$ unzip acflow.zip
```

or

```shell
$ tar xvfz acflow.tar.gz
```

Third, please input the following commands in Julia's REPL (Read-Eval-Print Loop) environment:

```julia-repl
julia> using Pkg
julia> Pkg.add(url = "/home/your_home/acflow/")
```

Here, `Pkg` is Julia's built-in package manager, and `/home/your_home/acflow` is assumed to be the root directory of the ACFlow toolkit. In practice, the users can use the `Pkg` package to install the ACFlow toolkit from its github repository directly:

```julia-repl
julia> using Pkg
julia> Pkg.add(url = "https://github.com/huangli712/ACFlow")
```

So the second step is optional. Furthermore, if the installed ACFlow toolkit is outdated, the users can use the following commands to upgrade ACFlow:

```julia
julia> using Pkg
julia> Pkg.update("ACFlow")
```

!!! note

    How to verify whether the installation is successful? Please input the following codes in Julia's REPL:

    ```julia-repl
    julia> using ACFlow
    ACFlow > Loading...
    ACFlow > Well, ACFlow is compiled and loaded (214 functions).
    ACFlow > We are ready to go!


    julia> welcome()
    ╔═╗╔═╗╔═╗┬  ┌─┐┬ ┬
    ╠═╣║  ╠╣ │  │ ││││
    ╩ ╩╚═╝╚  ┴─┘└─┘└┴┘
    A Modern Toolkit for Analytical Continuation Problems
    Package: ACFlow
    Version: 2.1.4-devel.250405
    Release: 2025/04


    ```

    If the `welcome()` function works correctly, it means that the ACFlow package is installed successfully.

Finally, in order to generate the documentation, please type the following commands in the terminal:

```shell
$ pwd
/home/your_home/acflow
$ cd docs
$ julia make.jl
```

After a few seconds, the documentation is built and saved in the `acflow/docs/build` directory if everything is OK. The home page of the documentation is `acflow/docs/build/index.html`. The users can open it with any web browsers.

## Technical Notes

!!! info

    Sometimes in offline environment, the package manager is disable. At this time, we should setup the environment variable `ACFLOW_HOME` as follows:

    ```shell
    $ export ACFLOW_HOME=/home/your_home/acflow/src
    ```

    And then plug the following code in front of your Julia's scripts:

    ```julia
    haskey(ENV,"ACFLOW_HOME") && pushfirst!(LOAD_PATH, ENV["ACFLOW_HOME"])
    ```

    Or else Julia won't find the ACFlow package.

!!! note

    After v1.8.0, the ACFlow toolkit relies on the `Zygote.jl` package to calculate gradient. The Zygota.jl package adopts the automatic differentiation approach, which is more efficient than the finite difference method. If the users have trouble in installing the Zygote.jl package, we provide an internal function to bypass it.

    Firstly, please disable Zygote in `ACFlow.jl` (about line 73).
    ```julia
    # using Zygote
    ```

    Secondly, fix `𝐽!(J::Vector{C64}, x::Vector{C64})` in `nac.jl` (about line 844).
    ```julia
    function 𝐽!(J::Vector{C64}, x::Vector{C64})
        #J .= Zygote.gradient(𝑓, x)[1]

        # Finite difference algorithm
        J .= gradient_via_fd(𝑓, x)
    end
    ```

    Note that `gradient_via_fd()` is based on the finite difference method, which is much slower and less accurate than the automatic differentiation approach. It is implemented in the `math.jl`.

!!! note

    If Julia fails to load the ACFlow package, please make sure: (1) it is correctly installed via `Pkg.add()`. (2) or the `ACFLOW_HOME` environment variable is correctly set and is put into `LOAD_PATH`.
