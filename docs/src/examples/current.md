# Current-Current Correlation Function

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
