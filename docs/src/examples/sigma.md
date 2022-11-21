

# Matsubara Self-Energy Function

In order to demonstrate usefulness of the ACFlow toolkit, four examples are illustrated in this section. These examples cover typical application scenarios of the ACFlow toolkit, including analytical continuations of Matsubara self-energy function, Matsubara Green's function, imaginary time Green's function, and current-current correlation function within the script mode or standard mode. All of the necessary source codes and data files, which can be used to reproduce the results as shown in this section, are placed in the \texttt{/home/your\_home/acflow/test/T*} folders.  

\subsection{Matsubara self-energy function\label{subsec:sigma}}

Now let us consider the following single-band Hubbard model on a Bethe lattice at first:
\begin{equation}
H = -t \sum_{\langle ij \rangle \sigma} c^{\dagger}_{i\sigma}c_{j\sigma}
 - \mu \sum_i n_i + U \sum_i n_{i\uparrow} n_{i\downarrow},
\end{equation}
where $t$ is the hopping parameter, $\mu$ is the chemical potential, $U$ is the Coulomb interaction, $n$ is the occupation number, $\sigma$ denotes the spin, $i$ and $j$ are site indices. This model is solved by using the dynamical mean-field theory (dubbed DMFT)~\cite{RevModPhys.68.13} with the hybridization expansion continuous-time quantum Monte Carlo solver (dubbed CT-HYB)~\cite{RevModPhys.83.349} as implemented in the $i$QIST package~\cite{HUANG2015140,HUANG2017423}. The parameters used in the DMFT + CT-HYB calculation are $t = 0.5$, $U = 2.0$, $\mu = 1.0$, and $\beta = 10.0$. Once the DMFT self-consistent calculation is finished, the Matsubara self-energy function $\Sigma(i\omega_n)$ is obtained. We are going to convert it to real frequency self-energy function $\Sigma(\omega)$. The data of Matsubara self-energy function $\Sigma(i\omega_n)$ have been preprocessed and stored in \texttt{siw.data}. This file contains five columns, which are used to record the Matsubara frequency $\omega_n$, Re$\Sigma(i\omega_n)$, Im$\Sigma(i\omega_n)$, error bar of Re$\Sigma(i\omega_n)$, error bar of Im$\Sigma(i\omega_n)$, respectively. Only the first twenty Matsubara frequency points are kept, because the high-frequency data are somewhat noisy.

The purpose of this example is to demonstrate usage of the \texttt{MaxEnt} solver and the script mode of the ACFlow toolkit. Next we will explain the key steps in detail. As for the complete Julia script, please refer to \texttt{sigma.jl} and \texttt{gendata.jl} in the \texttt{/home/your\_home/acflow/test/T01/} folder.   

First, we have to load the essential Julia packages. Both the \texttt{DelimitedFiles} and \texttt{Printf} packages belong to Julia's standard library. They are used to read input data and write calculated results, respectively.  

\begin{lstlisting}[language=Julia,
basicstyle=\ttfamily\small,
backgroundcolor=\color{yellow!10},
commentstyle=\color{olive!10!green},
keywordstyle=\color{purple}]
#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using DelimitedFiles
using Printf
using ACFlow

welcome() # Print welcome message only
\end{lstlisting}

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

\subsection{Matsubara Green's function\label{subsec:green}}

The purpose of the second example is to treat the Matsubara Green's function by using the \texttt{StochOM} solver.     

At first, please consider the following spectral density with two gaussian peaks:
\begin{equation}
A(\omega) = 
A_1 \exp\left[\frac{-(\omega - \epsilon_1)^2}{2 \Gamma^2_1}\right] +
A_2 \exp\left[\frac{-(\omega - \epsilon_2)^2}{2 \Gamma^2_2}\right],
\end{equation}
with $A_1 = 1.0$, $A_2 = 0.3$, $\epsilon_1 = 0.5$, $\epsilon_2 = -2.5$, $\Gamma_1 = 0.2$, and $\Gamma_2 = 0.8$. Then the Matsubara Green's function $G(i\omega_n)$ is evaluated by using Eq.~(\ref{eq:kernel_w}) with $\beta = 10.0$. Random noises, built by formula $0.0001 r_1 \exp(i 2\pi r_2 )$ where $r_1$ and $r_2$ are random numbers in (0.0,1.0), are added to $G(i\omega_n)$. The error bar of $G(i\omega_n)$ is fixed to 1e-4. The generated data for $G(i\omega_n)$ are written in \texttt{giw.data}.  

Next, we are going to use the standard mode, such that a configure file (\texttt{ac.toml}) must be prepared. It is listed as follows. Since the \texttt{StochOM} solver is chosen, the \texttt{[BASE]} and \texttt{[StochOM]} blocks must be present. 

\begin{lstlisting}[language=TOML,
basicstyle=\ttfamily\small,
backgroundcolor=\color{yellow!10},
commentstyle=\color{olive!10!green},
keywordstyle=\color{purple}]
[BASE]
finput = "giw.data"
solver = "StochOM"
ktype  = "fermi"
mtype  = "flat"
grid   = "ffreq"
mesh   = "linear"
ngrid  = 10
nmesh  = 501
wmax   = 5.0
wmin   = -5.0
beta   = 10.0
offdiag = false

[StochOM]
ntry  = 100000
nstep = 1000
nbox  = 100
sbox  = 0.005
wbox  = 0.02
norm  = -1.0
\end{lstlisting}

\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{T_E2.pdf}
\caption{Analytical continuation of Matsubara Green's function by using the stochastic optimization method. (a) Simulated and exact spectral functions. (b) Reconstructed and synthetic Matsubara Green's functions. Only the imaginary parts are presented in this figure. \label{fig:giw}}
\end{figure}

Then we use the \texttt{acrun.jl} or \texttt{Pacrun.jl} script to perform analytical continuation simulation. The calculated results are shown in Fig.~\ref{fig:giw}. As is seen in Fig.~\ref{fig:giw}(a), both the sharp peak around 0.5~eV and the broad peak around -2.5~eV are correctly reproduced by the \texttt{StochOM} solver. In Fig.~\ref{fig:giw}(b), the reconstructed Matsubara Green's function agrees quite well with the raw input data.

\subsection{Imaginary time Green's function\label{subsec:gtime}}

In this example, analytical continuation of imaginary time Green's function will be tested. Note that this example is borrowed from Reference~\cite{beach} directly.

The exact spectral function reads:
\begin{equation}
A(\omega) =
\begin{cases}
\frac{1}{W} \frac{|\omega|}{\sqrt{\omega^2 - \Delta^2}},~\quad & \text{if}~\Delta < |\omega| < W/2. \\
0, & \text{otherwise}.
\end{cases} 
\end{equation}
Here, $W$ denotes bandwidth, and $\Delta$ is used to control size of the energy gap. Let $W = 6$ and $2\Delta = 1$. This spectrum should exhibit flat shoulders, steep peaks, and sharp gap edges. Actually, it is the spectrum of a BCS superconductor. 

First, the imaginary time Green's function $G(\tau)$ is generated using Eq.~(\ref{eq:kernel_t}). Then a normally-distributed random noise is add to $G(\tau)$. Maximum amplitude of the noise is 1e-4. The error bar of $G(\tau)$ is fixed to 1e-3. The data are written in \texttt{gtau.data}.

Next, we try to prepare the configure file (\texttt{ac.toml}). In this case, we would like to benchmark the \texttt{StochAC} solver, so the \texttt{solver} parameter is set to ``StochAC'' and the \texttt{grid} parameter is set to ``ftime''. Furthermore, the \texttt{exclude} parameter is enabled to impose some \emph{a priori} constraints to the spectrum. The full \texttt{ac.toml} is listed as follows:

\begin{lstlisting}[language=TOML,
basicstyle=\ttfamily\small,
backgroundcolor=\color{yellow!10},
commentstyle=\color{olive!10!green},
keywordstyle=\color{purple}]
[BASE]
finput = "giw.data"
solver = "MaxEnt"
ktype  = "fermi"
mtype  = "flat"
grid   = "ffreq"
mesh   = "linear"
ngrid  = 10
nmesh  = 501
wmax   = 5.0
wmin   = -5.0
beta   = 10.0
offdiag = false
exclude = [[-5.0,-3.0], [-0.5,0.5], [3.0,5.0]]

[StochAC]
nfine = 10000
ngamm = 512
nwarm = 4000
nstep = 10000000
ndump = 40000
nalph = 40
alpha = 1.00
ratio = 1.20
\end{lstlisting}

We perform analytical continuation simulation by running the \texttt{acrun.jl} or \texttt{Pacrun.jl} script. In order to obtain smooth spectral density, it is useful to increase number of $\delta$ functions (See \texttt{ngamm} parameter) and number of Monte Carlo sampling steps (See \texttt{nstep} parameter).     

\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{T_E3.pdf}
\caption{Analytical continuation of imaginary time Green's function by using the stochastic analytical continuation (Beach's algorithm). (a) Simulated and exact spectral functions. (b) $\alpha$-dependent spectral functions. (c) Internal energy $U$ as a function of $\alpha$. The vertical bar indicates the optimal $\alpha$ parameter. (d) Simulated and exact imaginary time Green's functions. \label{fig:gtau}}
\end{figure}
 
Figure~\ref{fig:gtau} shows the calculated results. In Fig.~\ref{fig:gtau}(a), the exact spectral function is compared with the simulated spectrum. Note that besides the \texttt{StochAC} solver, the other three solvers are also tested. Their results are also plotted in this figure for a direct comparison. It is remarkable that the \texttt{StochAC} and \texttt{StochSK} solvers do a superior job of modelling the spectrum. The major characteristics of the spectrum, including flat regions, steep peaks, and sharp gap edges, are well captured by the two solvers. Especially, we have finished more tests without any constraints on the spectral density. The gap in the spectrum can be reproduced as well. On the other hand, the spectra obtained by the \texttt{MaxEnt} and \texttt{StochOM} solvers are much too smooth, and show extra shoulder peaks around $\pm$2.0. Figure~\ref{fig:gtau}(b) shows $\alpha$-resolved spectral functions $A_{\alpha}(\omega)$ for selected $\alpha$ parameters. Fluctuation in the flat regions of the calculated spectral density grows when $\alpha$ increases. Figure~\ref{fig:gtau}(c) shows internal energy $U$ as a function of $\alpha$. From this figure, the critical $\alpha$ is estimated, which is indicated by the vertical bar. Finally, the reproduced Green's function $\tilde{G}(\tau)$ agrees quite well with the raw input data, which is shown in Fig.~\ref{fig:gtau}(d).     

\subsection{Current-current correlation function\label{subsec:optic}}

The former three examples only concern fermionic correlators. How about bosonic correlators? In this example, we will demonstrate how to perform analytical continuation simulation for a typical bosonic correlator, the current-current correlation function $\Pi(\tau)$, to obtain the optical conductivity $\sigma(\omega)$. Note that this example is taken from Reference~\cite{PhysRevB.82.165125} directly. 

The exact optical conductivity $\sigma(\omega)$ reads:
\begin{equation}
\label{eq:optic}
\sigma(\omega) = 
\left\{
\frac{W_1}{1 + (\omega/\Gamma_1)^2} + 
\frac{W_2}{1 + [(\omega - \epsilon)/\Gamma_2]^2} +
\frac{W_2}{1 + [(\omega + \epsilon)/\Gamma_2]^2}
\right\}
\frac{1}{1 + (\omega/\Gamma_3)^6},
\end{equation}
where $W_1 = 0.3$, $W_2 = 0.2$, $\Gamma_1 = 0.3$, $\Gamma_2 = 1.2$, $\Gamma_3 = 4.0$, and $\epsilon = 3.0$. The current-current correlation function $\Pi(\tau)$ can be evaluated from $\sigma(\omega)$ by using the following equation:
\begin{equation}
\label{eq:current}
\Pi(\tau) = \int^{\infty}_{-\infty} K(\tau,\omega) \sigma(\omega)~d\omega,
\end{equation}
where the kernel function $K(\tau,\omega)$ is different from Eq.~(\ref{eq:ktau}). It reads:
\begin{equation}
\label{eq:Koptic}
K(\tau,\omega) = \frac{1}{\pi} \frac{\omega e^{-\tau\omega}}{1- e^{-\beta\omega}}.
\end{equation}
In this case, $\beta$ is fixed to be 20.0. 

At first, we use Eq.~(\ref{eq:optic}) $\sim$ Eq.~(\ref{eq:Koptic}) to prepare $\Pi(\tau)$. The error bar of $\Pi(\tau)$ is fixed to 1e-4. The calculated $\Pi(\tau)$ is written in \texttt{chit.data}. 

Next, we conduct analytical continuation simulation as usual. The used configuration file is attached as follows. Here, the \texttt{StochSK} solver is adopted, so the \texttt{solver} parameter is ``StochSK'' and the \texttt{grid} parameter is ``btime''. And the Shao-Sandvik algorithm~\cite{PhysRevX.7.041072} is applied to seek optimal $\Theta$, so the \texttt{method} parameter is ``chi2min''. The users can further increase the values of \texttt{nfine}, \texttt{ngamm}, and \texttt{nstep} parameters to improve computational accuracy. 
 
\begin{lstlisting}[language=TOML,
basicstyle=\ttfamily\small,
backgroundcolor=\color{yellow!10},
commentstyle=\color{olive!10!green},
keywordstyle=\color{purple}]
[BASE]
finput = "chit.data"
solver = "StochSK"
ktype  = "bsymm"
mtype  = "flat"
grid   = "btime"
mesh   = "halflorentz"
ngrid  = 10
nmesh  = 801
wmax   = 8.0
wmin   = 0.0
beta   = 20.0
offdiag = false

[StochSK]
method = "chi2min"
nfine = 40000
ngamm = 1000
nwarm = 1000
nstep = 20000
ndump = 200
retry = 10
theta = 1e+6
ratio = 0.90
\end{lstlisting}

\begin{figure}[ht]
\centering
\includegraphics[width=\textwidth]{T_E4.pdf}
\caption{Analytical continuation of current-current correlation function by using the stochastic analytical continuation (Sandvik's algorithm). (a) Simulated and exact optical conductivities $\sigma(\omega)$. (b) Simulated and exact current-current correlation functions $\Pi(\tau)$. \label{fig:optic}}
\end{figure}

The calculated results are illustrated in Fig.~\ref{fig:optic}. From Fig.~\ref{fig:optic}(a), it is clear that the main features of optical conductivity are successfully captured by the \texttt{StochSK} solver. Both the sharp Drude peak at $\omega = 0$ and a broad satellite peak around $\omega = 3.0$ are well reproduced. As is seen in Fig.~\ref{fig:optic}(b), the reconstructed $\tilde{\Pi}(\tau)$ coincides with the original $\Pi(\tau)$. 
