# Used for objectives and solvers where the gradient is available/exists
mutable struct OnceDifferentiable
    ℱ  # objective
    𝒥  # (partial) derivative of objective
    𝐹  # cache for f output
    𝐽  # cache for j output
end

function OnceDifferentiable(f, p0::AbstractArray, 𝐹::AbstractArray)
    function calc_F!(F, x)
        copyto!(F, f(x))
    end

    function calc_J!(J, x)
        relstep = cbrt(eps(real(eltype(x))))
        absstep = relstep
        @inbounds for i ∈ 1:length(x)
            x_save = x[i]
            epsilon = max(relstep * abs(x_save), absstep)
            x[i] = x_save + epsilon
            fx2 = vec(f(x))
            x[i] = x_save - epsilon
            fx1 = vec(f(x))
            J[:,i] = (fx2 - fx1) ./ (2 * epsilon)
            x[i] = x_save
        end
    end

    𝐽 = eltype(p0)(NaN) .* vec(𝐹) .* vec(p0)'
    OnceDifferentiable(calc_F!, calc_J!, 𝐹, 𝐽)
end

value(obj::OnceDifferentiable) = obj.𝐹
function value(obj::OnceDifferentiable, F, x)
    obj.ℱ(F, x)
end
function value!(obj::OnceDifferentiable, x)
    obj.ℱ(obj.𝐹, x)
    obj.𝐹
end

jacobian(obj::OnceDifferentiable) = obj.𝐽
function jacobian!(obj::OnceDifferentiable, x)
    obj.𝒥(obj.𝐽, x)
    obj.𝐽
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
    value!(df, initial_x)
    jacobian!(df, initial_x)
    
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
