# Used for objectives and solvers where the gradient is available/exists
mutable struct OnceDifferentiable
    ℱ  # objective
    𝒥  # (partial) derivative of objective
    𝐹  # cache for f output
    𝐽  # cache for j output
end

function OnceDifferentiable(𝑓, p0::AbstractArray, 𝐹::AbstractArray)
    function ℱ!(𝐹, x)
        copyto!(𝐹, 𝑓(x))
    end

    function 𝒥!(𝐽, x)
        rel_step = cbrt(eps(real(eltype(x))))
        abs_step = rel_step
        @inbounds for i ∈ 1:length(x)
            xₛ = x[i]
            ϵ = max(rel_step * abs(xₛ), abs_step)
            x[i] = xₛ + ϵ
            f₂ = vec(𝑓(x))
            x[i] = xₛ - ϵ
            f₁ = vec(𝑓(x))
            𝐽[:,i] = (f₂ - f₁) ./ (2 * ϵ)
            x[i] = xₛ
        end
    end

    𝐽 = eltype(p0)(NaN) .* vec(𝐹) .* vec(p0)'
    OnceDifferentiable(ℱ!, 𝒥!, 𝐹, 𝐽)
end

value(obj::OnceDifferentiable) = obj.𝐹
value(obj::OnceDifferentiable, 𝐹, x) = obj.ℱ(𝐹, x)
function value!(obj::OnceDifferentiable, x)
    obj.ℱ(obj.𝐹, x)
    obj.𝐹
end

jacobian(obj::OnceDifferentiable) = obj.𝐽
jacobian(obj::OnceDifferentiable, 𝐽, x) = obj.𝒥(𝐽, x)
function jacobian!(obj::OnceDifferentiable, x)
    obj.𝒥(obj.𝐽, x)
    obj.𝐽
end

mutable struct OptimizationResults{T,N}
    x₀::Array{T,N}
    minimizer::Array{T,N}
    minimum::T
    iterations::Int
    iteration_converged::Bool
    x_converged::Bool
    g_converged::Bool
end

"""
    levenberg_marquardt(df::OnceDifferentiable, x₀::AbstractVector{T})

Returns the argmin over x of `sum(f(x).^2)` using the Levenberg-Marquardt
algorithm, and an estimate of the Jacobian of `f` at x. The function `f`
is encoded in `df`. `x₀` is an initial guess for the solution.

See also: [`OnceDifferentiable`](@ref).
"""
function levenberg_marquardt(df::OnceDifferentiable, x₀::AbstractVector{T}) where T
    # Some constants
    max_lambda = 1e16 # minimum trust region radius
    min_lambda = 1e-16 # maximum trust region radius
    min_diagonal = 1e-6 # lower bound on values of diagonal matrix used to regularize the trust region step
    x_tol = 1e-8 # search tolerance in x
    g_tol = 1e-12 # search tolerance in gradient
    maxIter = 1000 # maximum number of iterations
    lambda = T(10) # (inverse of) initial trust region radius
    lambda_increase = 10.0 # lambda is multiplied by this factor after step below min quality
    lambda_decrease = 0.1 # lambda is multiplied by this factor after good quality steps
    min_step_quality = 1e-3 # for steps below this quality, the trust region is shrinked
    good_step_quality = 0.75 # for steps above this quality, the trust region is expanded

    # First evaluation
    value!(df, x₀)
    jacobian!(df, x₀)
    𝐹 = value(df)
    𝐽 = jacobian(df)

    converged = false
    x_converged = false
    g_converged = false
    iter = 0
    x = copy(x₀)

    trial_f = similar(𝐹)
    residual = sum(abs2, 𝐹)

    # Create buffers
    n = length(x)
    JJ = Matrix{T}(undef, n, n)
    n_buffer = Vector{T}(undef, n)
    Jdelta_buffer = similar(𝐹)

    while (~converged && iter < maxIter)
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

        DtD = vec(sum(abs2, 𝐽, dims=1))
        for i in 1:length(DtD)
            if DtD[i] <= min_diagonal
                DtD[i] = min_diagonal
            end
        end

        # delta_x = ( J'*J + lambda * Diagonal(DtD) ) \ ( -J'*F )
        mul!(JJ, 𝐽', 𝐽)
        @simd for i in 1:n
            @inbounds JJ[i, i] += lambda * DtD[i]
        end
        mul!(n_buffer, 𝐽', 𝐹)
        rmul!(n_buffer, -1)
        delta_x = JJ \ n_buffer

        # if the linear assumption is valid, our new residual should be:
        mul!(Jdelta_buffer, 𝐽, delta_x)
        Jdelta_buffer .= Jdelta_buffer .+ 𝐹
        predicted_residual = sum(abs2, Jdelta_buffer)

        # try the step and compute its quality
        # re-use n_buffer
        @. n_buffer = x + delta_x
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
            value!(df, n_buffer)
            residual = trial_residual
            if rho > good_step_quality
                # increase trust region radius
                lambda = max(lambda_decrease*lambda, min_lambda)
            end
        else
            # decrease trust region radius
            lambda = min(lambda_increase*lambda, max_lambda)
        end

        iter += 1

        # check convergence criteria:
        # 1. Small gradient: norm(J^T * value(df), Inf) < g_tol
        # 2. Small step size: norm(delta_x) < x_tol
        if norm(𝐽' * 𝐹, Inf) < g_tol
            g_converged = true
        end
        if norm(delta_x) < x_tol*(x_tol + norm(x))
            x_converged = true
        end
        converged = g_converged | x_converged
    end

    OptimizationResults(
        x₀,             # x₀
        x,                     # minimizer
        sum(abs2, value(df)),  # minimum
        iter,                # iterations
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
