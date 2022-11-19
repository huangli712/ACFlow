# Implementations

The ACFlow toolkit is developed with pure Julia language. Thanks to powerful type system and multiple dispatch paradigm of the Julia language, the four different analytical continuation solvers are integrated into an united software architecture. Redundant codes are greatly reduced. It is quite easy to implement new analytical continuation solver or add new features to the existing solvers in the future. Distributed computing is a built-in feature of Julia. So, it is straightforward to realize parallel calculations in the ACFlow toolkit. Now except for the \texttt{MaxEnt} solver, all the other solvers are paralleled.

\begin{table}[ht]
\centering
\begin{tabular}{l|l}
\hline
Filename & Description \\
\hline
\texttt{ACFlow.jl} & Entry of the ACFlow module. \\
\texttt{maxent.jl} & Maxent entropy method. \\
\texttt{sac.jl}    & Stochastic analytical continuation (K. S. D. Beach's version). \\
\texttt{san.jl}    & Stochastic analytical continuation (A. W. Sandvik's version). \\
\texttt{som.jl}    & Stochastic optimization method. \\
\texttt{global.jl} & Numerical and physical constants. \\
\texttt{types.jl}  & Basic data structures and computational parameters. \\
\texttt{base.jl}   & Driver for analytical continuation simulation. \\
\texttt{inout.jl}  & Read input data and write calculated results. \\
\texttt{config.jl} & Parse configuration file and extract computational parameters. \\
\texttt{math.jl}   & Root finding, numerical integration, interpolation, Einstein summation, and curve fitting. \\
\texttt{util.jl}   & Some utility functions. \\
\texttt{mesh.jl}   & Meshes for spectral density. \\
\texttt{grid.jl}   & Grids for input data. \\
\texttt{model.jl}  & Default model functions. \\
\texttt{kernel.jl} & Kernel functions. \\
\hline
\end{tabular}
\caption{List of source codes of the ACFlow toolkit. \label{tab:source}}
\end{table}

The source codes of the ACFlow toolkit are placed in the \texttt{acflow/src} folder. Their functions are summarized in Table~\ref{tab:source}. The documentation of the ACFlow toolkit is developed by using the Markdown language and \texttt{Documenter.jl} package. The source codes are placed in the \texttt{acflow/docs} folder. The users can build documentation by themselves. Please see Section~\ref{sec:usage} for how to do that. Or they can read the latest documentation in the following website:
\begin{verbatim}
    http://huangli712.github.io
\end{verbatim}    
Ten tests and four tutorials are also shipped with the ACFlow toolkit. Their source codes are placed in the \texttt{acflow/test} folder. See \texttt{acflow/test/test.md} and \texttt{acflow/test/tutor.md} for more details. The code repository of the ACFlow toolkit is:
\begin{verbatim}
    https://github.com/huangli712/ACFlow
\end{verbatim}