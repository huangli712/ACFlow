# Installation

It is an easy task to install the ACFlow toolkit. First, since it is written in pure Julia language, it is necessary to install the Julia runtime environment at first. The newest version of Julia is always preferred (version $>$ 1.60). Since the core codes only rely on Julia's built-in standard library, no the third-party packages are needed. Second, just download source codes of the ACFlow toolkit from its github repository. It should be a compressed file, such as `acflow.zip` or `acflow.tar.gz`. Please uncompress it in your favorite directory by using the following commands:

>    $ unzip acflow.zip

or

>    $ tar xvfz acflow.tar.gz

Third, the users have to announce a new environment variable `ACFLOW_HOME` for operation system. Supposed that the root directory of the ACFLow toolkit is `/home/your_home/acflow`, then `ACFLOW_HOME` should be setup as follows:

>    $ export ACFLOW_HOME=/home/your_home/acflow/src

Finally, in order to generate the documentation, the users should type the following commands in the terminal:

>    $ pwd
>    /home/your_home/acflow
>    $ cd docs
>    $ julia make.jl

After a few seconds, the documentation is built and saved in the `acflow/docs/build` directory if everything is OK. The home page of the documentation is `acflow/docs/build/index.html`. We can read it with any web browsers.