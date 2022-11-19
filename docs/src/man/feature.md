# Main Features

Now the ACFlow toolkit supports three analytical continuation methods as introduced above. It includes four different analytical continuation solvers, namely \texttt{MaxEnt}, \texttt{StochAC}, \texttt{StochSK}, and \texttt{StochOM}. Just as their names suggest, the \texttt{MaxEnt} solver implements the maximum entropy method~\cite{JARRELL1996133}. The \texttt{StochAC} and \texttt{StochSK} solvers implement the K. S. D. Beach's version~\cite{beach} and A. W. Sandvik's version~\cite{PhysRevB.57.10287} of stochastic analytical continuation, respectively. The \texttt{StochOM} solver implements the stochastic optimization method~\cite{PhysRevB.62.6317}. The ACFlow toolkit also provides a convenient library, which can be used to prepare and carry out analytical continuation calculations flexibly. The major features of the ACFlow toolkit are summarized in Table~\ref{tab:feature}.

\begin{table}[ht]
\centering
\begin{tabular}{l|l|l|l|l}
\hline
Features & \texttt{MaxEnt} & \texttt{StochAC} & \texttt{StochSK} & \texttt{StochOM} \\
\hline
Matrix-valued Green's function & Y & N & N & N \\
Imaginary time grid            & Y & Y & Y & Y \\
Matsubara frequency grid       & Y & Y & Y & Y \\
Linear mesh                    & Y & Y & Y & Y \\
Nonlinear mesh                 & Y & Y & Y & Y \\
Fermionic kernel               & Y & Y & Y & Y \\
Bosonic kernel                 & Y & Y & Y & Y \\
Self-defined model function    & Y & N & N & N \\
Constrained analytical continuation & N & Y & Y & Y \\
Regeneration of input data     & Y & Y & Y & Y \\
Kramers-Kronig transformation  & Y & Y & Y & Y \\
Parallel computing             & N & Y & Y & Y \\
Parallel tempering             & N & Y & N & N \\
Interactive mode               & Y & Y & Y & Y \\
Script mode                    & Y & Y & Y & Y \\
Standard mode                  & Y & Y & Y & Y \\
\hline
\end{tabular}
\caption{Major features of the ACFlow toolkit. \texttt{MaxEnt}, \texttt{StochAC}, \texttt{StochSK}, and \texttt{StochOM} are the four analytical continuation solvers implemented in this toolkit. \label{tab:feature}}
\end{table}

In Table~\ref{tab:feature}, ``Y'' means yes while ``N'' means no. ``Interactive mode'', ``Script mode'', and ``Standard model'' are three running modes supported by the ACFlow toolkit. We will introduce them in next section. The \texttt{MaxEnt} solver supports the ``historic'', ``classic'', ``bryan'', and ``chi2kink'' algorithms to determine the $\alpha$ parameter. The \texttt{StochAC} solver is only compatible with a flat model function, while the \texttt{StochSK} and \texttt{StochOM} solvers don't rely on any default model functions. The \texttt{StochOM} solver does not support analytical continuation of fermionic imaginary time Green's function for the moment. 