abstract type AbstractObjective end

# Initialize an n-by-n Jacobian
alloc_DF(x, F) = eltype(x)(NaN) .* vec(F) .* vec(x)'

x_of_nans(x, Tf=eltype(x)) = fill!(Tf.(x), Tf(NaN))

function f!_from_f(f, F::AbstractArray, inplace)
    if inplace
        return f
    else
        return function ff!(F, x)
            copyto!(F, f(x))
            F
        end
    end
end
function df!_from_df(j, F::AbstractArray, inplace)
    if inplace
        return j
    else
        return function jj!(J, x)
            copyto!(J, j(x))
            J
        end
    end
end
function fdf!_from_fdf(fj, F::AbstractArray, inplace)
    if inplace
        return fj
    else
        return function ffjj!(F, J, x)
            fx, jx = fj(x)
            copyto!(J, jx)
            copyto!(F, fx)
        end
    end
end

# Used for objectives and solvers where the gradient is available/exists
mutable struct OnceDifferentiable{TF, TDF, TX} <: AbstractObjective
    f # objective
    df # (partial) derivative of objective
    fdf # objective and (partial) derivative of objective
    F::TF # cache for f output
    DF::TDF # cache for df output
    x_f::TX # x used to evaluate f (stored in F)
    x_df::TX # x used to evaluate df (stored in DF)
    f_calls::Vector{Int}
    df_calls::Vector{Int}
end

### Only the objective
function OnceDifferentiable(f, x::AbstractArray, F::AbstractArray, DF::AbstractArray = alloc_DF(x, F); inplace = true, autodiff = :finite)
    f! = f!_from_f(f, F, inplace)
    OnceDifferentiable(f!, x::AbstractArray, F::AbstractArray, DF, autodiff)
end

function OnceDifferentiable(f, x_seed::AbstractArray, F::AbstractArray, DF::AbstractArray, autodiff::Symbol)
    # Figure out which Val-type to use for FiniteDiff based on our symbol interface.
    fdtype = Val{:central}
    # Apparently only the third input is aliased.
    j_finitediff_cache = FiniteDiff.JacobianCache(copy(x_seed), copy(F), copy(F), fdtype)

    function fj_finitediff!(F, J, x)
        f(F, x)
        FiniteDiff.finite_difference_jacobian!(J, f, x, j_finitediff_cache)
        F
    end
    function j_finitediff!(J, x)
        F_cache = copy(F)
        fj_finitediff!(F_cache, J, x)
    end

    return OnceDifferentiable(f, j_finitediff!, fj_finitediff!, x_seed, F, DF)
end

### Objective and derivative
function OnceDifferentiable(f, df, fdf, x::AbstractArray, F::AbstractArray, DF::AbstractArray = alloc_DF(x, F); inplace = true)
    f = f!_from_f(f, F, inplace)
    df! = df!_from_df(df, F, inplace)
    fdf! = fdf!_from_fdf(fdf, F, inplace)
    x_f, x_df = x_of_nans(x), x_of_nans(x)
    OnceDifferentiable(f, df!, fdf!, copy(F), copy(DF), x_f, x_df, [0,], [0,])
end

abstract type Optimizer end
struct LevenbergMarquardt <: Optimizer end
struct Avv
        h!::Function
    hessians::Array{Float64}

    function Avv(h!::Function, n::Int, m::Int)
        hessians = Array{Float64}(undef, m*n, n)
        new(h!, hessians)
    end
end
value(obj::AbstractObjective) = obj.F
function value(obj::AbstractObjective, F, x)
    obj.f_calls .+= 1
    return obj.f(F, x)
end
function value!(obj::AbstractObjective, x)
    if x != obj.x_f
        value!!(obj, x)
    end
    value(obj)
end
jacobian(obj::AbstractObjective) = obj.DF
value_jacobian!!(obj, x) = value_jacobian!!(obj, obj.F, obj.DF, x)
function value_jacobian!!(obj, F, J, x)
    obj.fdf(F, J, x)
    copyto!(obj.x_f, x)
    copyto!(obj.x_df, x)
    obj.f_calls .+= 1
    obj.df_calls .+= 1
    obj.df_calls
    F, J
end


function jacobian!(obj, x)
    if x != obj.x_df
        jacobian!!(obj, x)
    end
    jacobian(obj)
end

jacobian!!(obj, x) = jacobian!!(obj, obj.DF, x)
function jacobian!!(obj, J, x)
    obj.df(J, x)
    copyto!(obj.x_df, x)
    obj.df_calls .+= 1
    obj.df_calls
    J
end

struct OptimizationState{T <: Optimizer}
    iteration::Int
    value::Float64
    g_norm::Float64
    metadata::Dict
end
#=
function Base.show(io::IO, t::OptimizationState)
    @printf io "%6d   %14e   %14e\n" t.iteration t.value t.g_norm
    if !isempty(t.metadata)
        for (key, value) in t.metadata
            @printf io " * %s: %s\n" key value
        end
    end
    return
end
=#

OptimizationTrace{T} = Vector{OptimizationState{T}}
#=
function Base.show(io::IO, tr::OptimizationTrace)
    @printf io "Iter     Function value   Gradient norm \n"
    @printf io "------   --------------   --------------\n"
    for state in tr
        show(io, state)
    end
    return
end
=#

abstract type OptimizationResults end
minimizer(r::OptimizationResults) = r.minimizer
mutable struct MultivariateOptimizationResults{O<:Optimizer,T,N} <: OptimizationResults
    method::O
    initial_x::Array{T,N}
    minimizer::Array{T,N}
    minimum::T
    iterations::Int
    iteration_converged::Bool
    x_converged::Bool
    x_tol::Float64
    x_residual::Float64
    f_converged::Bool
    f_tol::Float64
    f_residual::Float64
    g_converged::Bool
    g_tol::Float64
    g_residual::Float64
    f_increased::Bool
    trace::OptimizationTrace{O}
    f_calls::Int
    g_calls::Int
    h_calls::Int
end
converged(r::MultivariateOptimizationResults) = r.x_converged || r.f_converged || r.g_converged

Base.summary(::LevenbergMarquardt) = "Levenberg-Marquardt"
"""
    `levenberg_marquardt(f, g, initial_x; <keyword arguments>`

Returns the argmin over x of `sum(f(x).^2)` using the Levenberg-Marquardt
algorithm, and an estimate of the Jacobian of `f` at x.

The function `f` should take an input vector of length n and return an output
vector of length m. The function `g` is the Jacobian of f, and should return an m x
n matrix. `initial_x` is an initial guess for the solution.

Implements box constraints as described in Kanzow, Yamashita, Fukushima (2004; J
Comp & Applied Math).

# Keyword arguments
* `x_tol::Real=1e-8`: search tolerance in x
* `g_tol::Real=1e-12`: search tolerance in gradient
* `maxIter::Integer=1000`: maximum number of iterations
* `min_step_quality=1e-3`: for steps below this quality, the trust region is shrinked
* `good_step_quality=0.75`: for steps above this quality, the trust region is expanded
* `lambda::Real=10`: (inverse of) initial trust region radius
* `tau=Inf`: set initial trust region radius using the heuristic : tau*maximum(jacobian(df)'*jacobian(df))
* `lambda_increase=10.0`: `lambda` is multiplied by this factor after step below min quality
* `lambda_decrease=0.1`: `lambda` is multiplied by this factor after good quality steps
* `show_trace::Bool=false`: print a status summary on each iteration if true
* `lower,upper=[]`: bound solution to these limits
"""

# I think a smarter way to do this *might* be to create a type similar to `OnceDifferentiable`
# and the like. This way we could not only merge the two functions, but also have a convinient
# way to provide an autodiff-made acceleration when someone doesn't provide an `avv`.
# it would probably be very inefficient performace-wise for most cases, but it wouldn't hurt to have it somewhere
function levenberg_marquardt(df::OnceDifferentiable, initial_x::AbstractVector{T};
    x_tol::Real = 1e-8, g_tol::Real = 1e-12, maxIter::Integer = 1000,
    lambda = T(10), tau=T(Inf), lambda_increase::Real = 10.0, lambda_decrease::Real = 0.1,
    min_step_quality::Real = 1e-3, good_step_quality::Real = 0.75,
    show_trace::Bool = false, lower::Vector{T} = Array{T}(undef, 0), upper::Vector{T} = Array{T}(undef, 0), avv!::Union{Function,Nothing,Avv} = nothing
    ) where T

    # First evaluation
    value_jacobian!!(df, initial_x)
    
    if isfinite(tau)
        lambda = tau*maximum(jacobian(df)'*jacobian(df))
    end


    # check parameters
    ((isempty(lower) || length(lower)==length(initial_x)) && (isempty(upper) || length(upper)==length(initial_x))) ||
            throw(ArgumentError("Bounds must either be empty or of the same length as the number of parameters."))
    ((isempty(lower) || all(initial_x .>= lower)) && (isempty(upper) || all(initial_x .<= upper))) ||
            throw(ArgumentError("Initial guess must be within bounds."))
    (0 <= min_step_quality < 1) || throw(ArgumentError(" 0 <= min_step_quality < 1 must hold."))
    (0 < good_step_quality <= 1) || throw(ArgumentError(" 0 < good_step_quality <= 1 must hold."))
    (min_step_quality < good_step_quality) || throw(ArgumentError("min_step_quality < good_step_quality must hold."))


    # other constants
    MAX_LAMBDA = 1e16 # minimum trust region radius
    MIN_LAMBDA = 1e-16 # maximum trust region radius
    MIN_DIAGONAL = 1e-6 # lower bound on values of diagonal matrix used to regularize the trust region step


    converged = false
    x_converged = false
    g_converged = false
    iterCt = 0
    x = copy(initial_x)
    delta_x = copy(initial_x)
    a = similar(x)

    trial_f = similar(value(df))
    residual = sum(abs2, value(df))

    # Create buffers
    n = length(x)
    m = length(value(df))
    JJ = Matrix{T}(undef, n, n)
    n_buffer = Vector{T}(undef, n)
    Jdelta_buffer = similar(value(df))

    # and an alias for the jacobian
    J = jacobian(df)
    dir_deriv = Array{T}(undef,m)
    v = Array{T}(undef,n)

    # Maintain a trace of the system.
    tr = OptimizationTrace{LevenbergMarquardt}()
    if show_trace
        d = Dict("lambda" => lambda)
        os = OptimizationState{LevenbergMarquardt}(iterCt, sum(abs2, value(df)), NaN, d)
        push!(tr, os)
        println(os)
    end

    while (~converged && iterCt < maxIter)
        # jacobian! will check if x is new or not, so it is only actually
        # evaluated if x was updated last iteration.
        jacobian!(df, x) # has alias J

        # we want to solve:
        #    argmin 0.5*||J(x)*delta_x + f(x)||^2 + lambda*||diagm(J'*J)*delta_x||^2
        # Solving for the minimum gives:
        #    (J'*J + lambda*diagm(DtD)) * delta_x == -J' * f(x), where DtD = sum(abs2, J,1)
        # Where we have used the equivalence: diagm(J'*J) = diagm(sum(abs2, J,1))
        # It is additionally useful to bound the elements of DtD below to help
        # prevent "parameter evaporation".

        DtD = vec(sum(abs2, J, dims=1))
        for i in 1:length(DtD)
            if DtD[i] <= MIN_DIAGONAL
                DtD[i] = MIN_DIAGONAL
            end
        end

        # delta_x = ( J'*J + lambda * Diagonal(DtD) ) \ ( -J'*value(df) )
        mul!(JJ, transpose(J), J)
        @simd for i in 1:n
            @inbounds JJ[i, i] += lambda * DtD[i]
        end
        #n_buffer is delta C, JJ is g compared to Mark's code
        mul!(n_buffer, transpose(J), value(df))
        rmul!(n_buffer, -1)

        v .= JJ \ n_buffer


        if avv! != nothing
            #GEODESIC ACCELERATION PART
            avv!(dir_deriv, x, v)
            mul!(a, transpose(J), dir_deriv)
            rmul!(a, -1) #we multiply by -1 before the decomposition/division
            LAPACK.potrf!('U', JJ) #in place cholesky decomposition
            LAPACK.potrs!('U', JJ, a) #divides a by JJ, taking into account the fact that JJ is now the `U` cholesky decoposition of what it was before
            rmul!(a, 0.5)
            delta_x .= v .+ a
            #end of the GEODESIC ACCELERATION PART
        else
            delta_x = v
        end





        # apply box constraints
        if !isempty(lower)
            @simd for i in 1:n
               @inbounds delta_x[i] = max(x[i] + delta_x[i], lower[i]) - x[i]
            end
        end
        if !isempty(upper)
            @simd for i in 1:n
               @inbounds delta_x[i] = min(x[i] + delta_x[i], upper[i]) - x[i]
            end
        end

        # if the linear assumption is valid, our new residual should be:
        mul!(Jdelta_buffer, J, delta_x)
        Jdelta_buffer .= Jdelta_buffer .+ value(df)
        predicted_residual = sum(abs2, Jdelta_buffer)

        # try the step and compute its quality
        # compute it inplace according to NLSolversBase value(obj, cache, state)
        # interface. No bang (!) because it doesn't update df besides mutating
        # the number of f_calls

        # re-use n_buffer
        n_buffer .= x .+ delta_x
        value(df, trial_f, n_buffer)

        # update the sum of squares
        trial_residual = sum(abs2, trial_f)

        # step quality = residual change / predicted residual change
        rho = (trial_residual - residual) / (predicted_residual - residual)
        if rho > min_step_quality
            # apply the step to x - n_buffer is ready to be used by the delta_x
            # calculations after this step.
            x .= n_buffer
            # There should be an update_x_value to do this safely
            copyto!(df.x_f, x)
            copyto!(value(df), trial_f)
            residual = trial_residual
            if rho > good_step_quality
                # increase trust region radius
                lambda = max(lambda_decrease*lambda, MIN_LAMBDA)
            end
        else
            # decrease trust region radius
            lambda = min(lambda_increase*lambda, MAX_LAMBDA)
        end

        iterCt += 1

        # show state
        if show_trace
            g_norm = norm(J' * value(df), Inf)
            d = Dict("g(x)" => g_norm, "dx" => delta_x, "lambda" => lambda)
            os = OptimizationState{LevenbergMarquardt}(iterCt, sum(abs2, value(df)), g_norm, d)
            push!(tr, os)
            println(os)
        end

        # check convergence criteria:
        # 1. Small gradient: norm(J^T * value(df), Inf) < g_tol
        # 2. Small step size: norm(delta_x) < x_tol
        if norm(J' * value(df), Inf) < g_tol
            g_converged = true
        end
        if norm(delta_x) < x_tol*(x_tol + norm(x))
            x_converged = true
        end
        converged = g_converged | x_converged
    end

    MultivariateOptimizationResults(
        LevenbergMarquardt(),    # method
        initial_x,             # initial_x
        x,                     # minimizer
        sum(abs2, value(df)),       # minimum
        iterCt,                # iterations
        !converged,            # iteration_converged
        x_converged,           # x_converged
        0.0,                   # x_tol
        0.0,
        false,                 # f_converged
        0.0,                   # f_tol
        0.0,
        g_converged,           # g_converged
        g_tol,                  # g_tol
        0.0,
        false,                 # f_increased
        tr,                    # trace
        first(df.f_calls),               # f_calls
        first(df.df_calls),               # g_calls
        0                      # h_calls
    )
end

struct LsqFitResult{P, R, J, W <: AbstractArray}
    param::P
    resid::R
    jacobian::J
    converged::Bool
    wt::W
end

function check_data_health(xdata, ydata)
    if any(ismissing, xdata) || any(ismissing, ydata)
        error("Data contains `missing` values and a fit cannot be performed")
    end
    if any(isinf, xdata) || any(isinf, ydata) || any(isnan, xdata) || any(isnan, ydata)
        error("Data contains `Inf` or `NaN` values and a fit cannot be performed")
    end
end

function lmfit(f, p0::AbstractArray, wt::AbstractArray; autodiff = :finite, kwargs...)
    # this is a convenience function for the curve_fit() methods
    # which assume f(p) is the cost functionj i.e. the residual of a
    # model where
    #   model(xpts, params...) = ydata + error (noise)

    # this minimizes f(p) using a least squares sum of squared error:
    #   rss = sum(f(p)^2)
    #
    # returns p, f(p), g(p) where
    #   p    : best fit parameters
    #   f(p) : function evaluated at best fit p, (weighted) residuals
    #   g(p) : estimated Jacobian at p (Jacobian with respect to p)

    # construct Jacobian function, which uses finite difference method
    r = f(p0)
    autodiff = autodiff == :forwarddiff ? :forward : autodiff
    R = OnceDifferentiable(f, p0, copy(r); inplace = false, autodiff = autodiff)
    results = levenberg_marquardt(R, p0; kwargs...)
    p = minimizer(results)
    return LsqFitResult(p, value!(R, p), jacobian!(R, p), converged(results), wt)
end

"""
    curve_fit(model, xdata, ydata, p0) -> fit
Fit data to a non-linear `model`. `p0` is an initial model parameter guess (see Example).
The return object is a composite type (`LsqFitResult`), with some interesting values:

* `fit.resid` : residuals = vector of residuals
* `fit.jacobian` : estimated Jacobian at solution

additionally, it is possible to quiry the degrees of freedom with

* `dof(fit)`
* `coef(fit)`

## Example
```julia
# a two-parameter exponential model
# x: array of independent variables
# p: array of model parameters
model(x, p) = p[1]*exp.(-x.*p[2])

# some example data
# xdata: independent variables
# ydata: dependent variable
xdata = range(0, stop=10, length=20)
ydata = model(xdata, [1.0 2.0]) + 0.01*randn(length(xdata))
p0 = [0.5, 0.5]

fit = curve_fit(model, xdata, ydata, p0)
```
"""
function curve_fit(model, xdata::AbstractArray, ydata::AbstractArray, p0::AbstractArray; kwargs...)
    check_data_health(xdata, ydata)
    # construct the cost function
    T = eltype(ydata)

    f = (p) -> model(xdata, p) - ydata
    lmfit(f, p0, T[]; kwargs...)
end
