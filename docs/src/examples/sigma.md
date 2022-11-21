!!! info

In order to demonstrate usefulness of the ACFlow toolkit, four examples are illustrated in this section. These examples cover typical application scenarios of the ACFlow toolkit, including analytical continuations of Matsubara self-energy function, Matsubara Green's function, imaginary time Green's function, and current-current correlation function within the script mode or standard mode. All of the necessary source codes and data files, which can be used to reproduce the results as shown in this section, are placed in the `/home/your_home/acflow/test/T*` folders. 

# Matsubara Self-Energy Function

Now let us consider the following single-band Hubbard model on a Bethe lattice at first:
```math
\begin{equation}
H = -t \sum_{\langle ij \rangle \sigma} c^{\dagger}_{i\sigma}c_{j\sigma}
 - \mu \sum_i n_i + U \sum_i n_{i\uparrow} n_{i\downarrow},
\end{equation}
```
where $t$ is the hopping parameter, $\mu$ is the chemical potential, $U$ is the Coulomb interaction, $n$ is the occupation number, $\sigma$ denotes the spin, $i$ and $j$ are site indices. This model is solved by using the dynamical mean-field theory (dubbed DMFT)~\cite{RevModPhys.68.13} with the hybridization expansion continuous-time quantum Monte Carlo solver (dubbed CT-HYB)~\cite{RevModPhys.83.349} as implemented in the $i$QIST package~\cite{HUANG2015140,HUANG2017423}. The parameters used in the DMFT + CT-HYB calculation are $t = 0.5$, $U = 2.0$, $\mu = 1.0$, and $\beta = 10.0$. Once the DMFT self-consistent calculation is finished, the Matsubara self-energy function $\Sigma(i\omega_n)$ is obtained. We are going to convert it to real frequency self-energy function $\Sigma(\omega)$. The data of Matsubara self-energy function $\Sigma(i\omega_n)$ have been preprocessed and stored in \texttt{siw.data}. This file contains five columns, which are used to record the Matsubara frequency $\omega_n$, Re$\Sigma(i\omega_n)$, Im$\Sigma(i\omega_n)$, error bar of Re$\Sigma(i\omega_n)$, error bar of Im$\Sigma(i\omega_n)$, respectively. Only the first twenty Matsubara frequency points are kept, because the high-frequency data are somewhat noisy.

The purpose of this example is to demonstrate usage of the \texttt{MaxEnt} solver and the script mode of the ACFlow toolkit. Next we will explain the key steps in detail. As for the complete Julia script, please refer to \texttt{sigma.jl} and \texttt{gendata.jl} in the \texttt{/home/your\_home/acflow/test/T01/} folder.   

First, we have to load the essential Julia packages. Both the \texttt{DelimitedFiles} and \texttt{Printf} packages belong to Julia's standard library. They are used to read input data and write calculated results, respectively.  

>#!/usr/bin/env julia
>
>push!(LOAD_PATH, ENV["ACFLOW_HOME"])
>
>using DelimitedFiles
>
>using Printf
>
>using ACFlow
>
>welcome() # Print welcome message only

Next, the data of Matsubara self-energy function are read from \texttt{siw.data}. The Hartree term $\Sigma_{H}$ should be subtracted from its real part:
\begin{equation}
\Sigma(i\omega_n) \to \Sigma(i\omega_n) - \Sigma_{H}.
\end{equation}
Note that $\Sigma_{H}$ is approximately equal to the asymptotic value of real part of $\Sigma(i\omega_n)$ when $n$ goes to infinite.   
 
\begin{lstlisting}[language=Julia,
basicstyle=\ttfamily\small,
backgroundcolor=\color{yellow!10},
commentstyle=\color{olive!10!green},
keywordstyle=\color{purple}]
# Deal with self-energy function
#
# Read self-energy function
dlm = readdlm("siw.data")
#
# Get grid
grid = dlm[:,1]
#
# Get self-energy function
Sinp = dlm[:,2] + im * dlm[:,3] # Value
Serr = dlm[:,4] + im * dlm[:,5] # Error bar
#
# Subtract hartree term
Sh = 1.0
@. Sinp = Sinp - Sh
\end{lstlisting}

Next, the computational parameters are encapsulated into two dictionaries. Then the \texttt{setup\_param()} function is called, so that these parameters take effect. Here, the \texttt{MatEnt} solver~\cite{JARRELL1996133,PhysRevB.44.6011} is employed to tackle the analytical continuation problem. But the other stochastic sampling solvers are also applicable. The default model function is gaussian. The mesh for spectral density is non-uniform (A tangent mesh). The number of used $\alpha$ parameters is 15, and the optimal $\alpha$ parameter is determined by the $\chi^2$kink algorithm~\cite{PhysRevE.94.023303}. 

\begin{lstlisting}[language=Julia,
basicstyle=\ttfamily\small,
backgroundcolor=\color{yellow!10},
commentstyle=\color{olive!10!green},
keywordstyle=\color{purple}]
# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "solver" => "MaxEnt",  # Choose MaxEnt solver
    "mtype"  => "gauss",   # Default model function
    "mesh"   => "tangent", # Mesh for spectral density
    "ngrid"  => 20,        # Number of input points
    "nmesh"  => 801,       # Number of output points
    "wmax"   => 8.0,       # Right boundary of mesh
    "wmin"   => -8.0,      # Left boundary of mesh
    "beta"   => 10.0,      # Inverse temperature
)
#
# For [MaxEnt] block
# See types.jl/_PMaxEnt for default setup
S = Dict{String,Any}(
    "nalph"  => 15,        # Number of alpha
    "alpha"  => 1e12,      # Starting value of alpha
    "blur"   => -1.0,      # Enable preblur or not
)
#
# Let the parameters take effect
setup_param(B, S)
\end{lstlisting}

It is quite easy to start the analytical continuation calculation. Just call the \texttt{solve()} function and pass the grid, input data, and error bar data to it. The return values of this function call are real frequency mesh, spectral density, and reconstructed Matsubara self-energy function. 

\begin{lstlisting}[language=Julia,
basicstyle=\ttfamily\small,
backgroundcolor=\color{yellow!10},
commentstyle=\color{olive!10!green},
keywordstyle=\color{purple}]
# Call the solver
mesh, Aout, Sout = solve(grid, Sinp, Serr)
\end{lstlisting}

Finally, the real frequency self-energy function must be supplemented with the Hartree term. Then the final results are written into \texttt{sigma.data}.   
   
\begin{lstlisting}[language=Julia,
basicstyle=\ttfamily\small,
backgroundcolor=\color{yellow!10},
commentstyle=\color{olive!10!green},
keywordstyle=\color{purple}]
# Calculate final self-energy function on real axis
#
# Add hartree term
@. Sout = Sout + Sh
#
# Write self-energy function to sigma.data
open("sigma.data", "w") do fout
    for i in eachindex(mesh)
        z = Sout[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n",
                mesh[i], real(z), imag(z))
    end
end
\end{lstlisting}

\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{T_E1.pdf}
\caption{Analytical continuation of Matsubara self-energy function by using the maximum entropy method. (a) Real part of real frequency self-energy function. (b) Imaginary part of real frequency self-energy function. (c) $\chi^{2}$ as a function of $\alpha$. The vertical bar indicates the optimal $\alpha$ parameter chosen by the $\chi^2$kink algorithm. (d) Reproduced and original data for imaginary part of the Matsubara self-energy functions. \label{fig:sig}}
\end{figure}

The calculated results are displayed in Fig.~\ref{fig:sig}. Fig.~{\ref{fig:sig}}(a) and (b) show the real and imaginary parts of the real frequency self-energy function. Near the Fermi level, Re$\Sigma(\omega)$ exhibits quasi-linear behavior, with which the quasiparticle weight $Z$ and effective mass of electron $m^*$ can be easily evaluated. As for the imaginary part, Im$\Sigma(0)$ is finite, which indicates that the electron-electron scattering is not trivial. Fig.~\ref{fig:sig}(c) shows the $\alpha$-dependent $\chi^{2}$. The vertical bar in this figure indicates the optimal $\alpha$ is around 2.15. In Fig.~\ref{fig:sig}(d), the reproduced and raw self-energy functions are compared. It is apparent that they are consistent with each other.
