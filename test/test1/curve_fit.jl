# Used for objectives and solvers where the gradient is available/exists
mutable struct OnceDifferentiable
    f    # objective
    df   # (partial) derivative of objective
    fdf  # objective and (partial) derivative of objective
    F    # cache for f output
    DF   # cache for df output
    x_f  # x used to evaluate f (stored in F)
    x_df # x used to evaluate df (stored in DF)
end

function OnceDifferentiable(f, p0::AbstractArray, F::AbstractArray)
    function calc_jacobian!(J, f!, x)
        relstep = cbrt(eps(real(eltype(x))))
        absstep = relstep
        x1, fx, fx1 = cache
        copyto!(x1, x)
        vfx = vec(fx)
        vfx1 = vec(fx1)
        @inbounds for i ∈ 1:length(x1)
            x_save = x[i]
            epsilon = max(relstep * abs(x_save), absstep)
            x1[i] = x_save + epsilon
            f!(fx1, x1)
            x1[i] = x_save - epsilon
            f!(fx, x1)
            @. J[:,i] = (vfx1 - vfx) / (2 * epsilon)
            x1[i] = x_save
        end
        nothing
    end

    function calc_F!(F, x)
        copyto!(F, f(x))
    end

    function calc_J!(J, x)
        calc_jacobian!(J, calc_F!, x)
    end

    function calc_FJ!(F, J, x)
        copyto!(F, f(x))
        calc_jacobian!(J, calc_F!, x)
    end

    DF = eltype(p0)(NaN) .* vec(F) .* vec(p0)'
    cache = (copy(p0), copy(F), copy(F))
    x_f = fill(NaN, length(p0))
    x_df = fill(NaN, length(p0))
    OnceDifferentiable(calc_F!, calc_J!, calc_FJ!, F, DF, x_f, x_df)
end

value(obj::OnceDifferentiable) = obj.F
function value(obj::OnceDifferentiable, F, x)
    return obj.f(F, x)
end
function value!(obj::OnceDifferentiable, x)
    if x != obj.x_f
        value!!(obj, x)
    end
    value(obj)
end

jacobian(obj::OnceDifferentiable) = obj.DF
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
    J
end

value_jacobian!!(obj, x) = value_jacobian!!(obj, obj.F, obj.DF, x)
function value_jacobian!!(obj, F, J, x)
    obj.fdf(F, J, x)
    copyto!(obj.x_f, x)
    copyto!(obj.x_df, x)
    F, J
end

mutable struct OptimizationResults{T,N}
    initial_x::Array{T,N}
    minimizer::Array{T,N}
    minimum::T
    iterations::Int
    iteration_converged::Bool
    x_converged::Bool
    g_converged::Bool
end

"""
    levenberg_marquardt(f, initial_x; kwargs...)

Returns the argmin over x of `sum(f(x).^2)` using the Levenberg-Marquardt
algorithm, and an estimate of the Jacobian of `f` at x. The function `f`
should take an input vector of length n and return an output vector of
length m. `initial_x` is an initial guess for the solution.

* x_tol, search tolerance in x
* g_tol, search tolerance in gradient
* maxIter, maximum number of iterations
* lambda, (inverse of) initial trust region radius
* lambda_increase, lambda is multiplied by this factor after step below min quality
* lambda_decrease, lambda is multiplied by this factor after good quality steps
* min_step_quality, for steps below this quality, the trust region is shrinked
* good_step_quality, for steps above this quality, the trust region is expanded
"""
function levenberg_marquardt(df::OnceDifferentiable, initial_x::AbstractVector{T};
    x_tol::Real = 1e-8,
    g_tol::Real = 1e-12,
    maxIter::Integer = 1000,
    lambda = T(10),
    lambda_increase::Real = 10.0,
    lambda_decrease::Real = 0.1,
    min_step_quality::Real = 1e-3,
    good_step_quality::Real = 0.75
) where T

    # First evaluation
    value_jacobian!!(df, initial_x)
    
    # other constants
    max_lambda = 1e16 # minimum trust region radius
    min_lambda = 1e-16 # maximum trust region radius
    min_diagonal = 1e-6 # lower bound on values of diagonal matrix used to regularize the trust region step

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
    v = Array{T}(undef,n)

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
            if DtD[i] <= min_diagonal
                DtD[i] = min_diagonal
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

        delta_x = v

        # if the linear assumption is valid, our new residual should be:
        mul!(Jdelta_buffer, J, delta_x)
        Jdelta_buffer .= Jdelta_buffer .+ value(df)
        predicted_residual = sum(abs2, Jdelta_buffer)

        # try the step and compute its quality
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
                lambda = max(lambda_decrease*lambda, min_lambda)
            end
        else
            # decrease trust region radius
            lambda = min(lambda_increase*lambda, max_lambda)
        end

        iterCt += 1

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

    OptimizationResults(
        initial_x,             # initial_x
        x,                     # minimizer
        sum(abs2, value(df)),  # minimum
        iterCt,                # iterations
        !converged,            # iteration_converged
        x_converged,           # x_converged
        g_converged,           # g_converged
    )
end

struct LsqFitResult
    param
    resid
    jacobian
    converged
end

"""
    curve_fit(model, x, y, p0)

Fit data to a non-linear `model`. `p0` is an initial model parameter guess.
The return object is a composite type (`LsqFitResult`).
"""
function curve_fit(model, x::AbstractArray, y::AbstractArray, p0::AbstractArray)
    f = (p) -> model(x, p) - y
    r = f(p0)
    R = OnceDifferentiable(f, p0, r)
    results = levenberg_marquardt(R, p0)
    p = results.minimizer
    conv = results.x_converged || results.g_converged
    return LsqFitResult(p, value!(R, p), jacobian!(R, p), conv)
end
