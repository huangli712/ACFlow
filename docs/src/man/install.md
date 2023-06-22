# [Installation](@id install)

It is an easy task to install the ACFlow toolkit. First, since it is written in pure Julia language, it is necessary to install the Julia runtime environment in advance. The newest version of Julia is always preferred (version > 1.60). Because the core codes only rely on Julia's built-in standard library, no the third-party packages are needed. Second, just download source codes of the ACFlow toolkit from its github repository. It should be a compressed file, such as `acflow.zip` or `acflow.tar.gz`. Please uncompress it in your favorite directory by using the following commands:

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

Here, `Pkg` is Julia's built-in package manager, and `/home/your\_home/acflow` is assumed to be the root directory of the ACFlow toolkit. In practice, the users can use the `Pkg` package to install the ACFlow toolkit from its github repository directly:

```julia-repl
julia> using Pkg
julia> Pkg.add(url = "https://github.com/huangli712/ACFlow")
```

So the second step is optional. Furthermore, if the installed ACFlow toolkit is outdated, the users can use the following commands to upgrade ACFlow:

```julia-repl
julia> using Pkg
julia> Pkg.update("ACFlow")
```

Finally, in order to generate the documentation, please type the following commands in the terminal:

```shell
$ pwd
/home/your_home/acflow
$ cd docs
$ julia make.jl
```

After a few seconds, the documentation is built and saved in the `acflow/docs/build` directory if everything is OK. The home page of the documentation is `acflow/docs/build/index.html`. The users can open it with any web browsers.
