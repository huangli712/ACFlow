# Used for objectives and solvers where the gradient is available/exists
mutable struct OnceDifferentiable
    ℱ! # objective
    𝒥! # (partial) derivative of objective
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
value(obj::OnceDifferentiable, 𝐹, x) = obj.ℱ!(𝐹, x)
function value!(obj::OnceDifferentiable, x)
    obj.ℱ!(obj.𝐹, x)
    obj.𝐹
end

jacobian(obj::OnceDifferentiable) = obj.𝐽
jacobian(obj::OnceDifferentiable, 𝐽, x) = obj.𝒥!(𝐽, x)
function jacobian!(obj::OnceDifferentiable, x)
    obj.𝒥!(obj.𝐽, x)
    obj.𝐽
end

struct OptimizationResults{T,N}
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
    Λₘ = 1e16 # minimum trust region radius
    λₘ = 1e-16 # maximum trust region radius
    min_diagonal = 1e-6 # lower bound on values of diagonal matrix
    x_tol = 1e-8 # search tolerance in x
    g_tol = 1e-12 # search tolerance in gradient
    maxIter = 1000 # maximum number of iterations
    λ = T(10) # (inverse of) initial trust region radius
    λᵢ = 10.0 # λ is multiplied by this factor after step below min quality
    λᵣ = 0.1 # λ is multiplied by this factor after good quality steps
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
    C_resid = sum(abs2, 𝐹)

    # Create buffers
    𝐽ᵀ𝐽 = diagm(x)
    𝐽δx = similar(𝐹)

    while (~converged && iter < maxIter)
        # Update jacobian 𝐽
        jacobian!(df, x)

        # Solve the equation: [𝐽ᵀ𝐽 + λ diag(𝐽ᵀ𝐽)] δ = 𝐽ᵀ𝐹
        mul!(𝐽ᵀ𝐽, 𝐽', 𝐽)
        𝐷ᵀ𝐷 = diag(𝐽ᵀ𝐽)
        replace!(x -> x ≤ min_diagonal ? min_diagonal : x, 𝐷ᵀ𝐷)
        @simd for i in eachindex(𝐷ᵀ𝐷)
            @inbounds 𝐽ᵀ𝐽[i,i] += λ * 𝐷ᵀ𝐷[i]
        end
        δx = - 𝐽ᵀ𝐽 \ (𝐽' * 𝐹)

        # If the linear assumption is valid, our new residual should be:
        mul!(𝐽δx, 𝐽, δx)
        𝐽δx .= 𝐽δx .+ 𝐹
        P_resid = sum(abs2, 𝐽δx)

        # Try to calculate new x, and then 𝐹, and then the residual.
        xnew = x + δx
        value(df, trial_f, xnew)
        T_resid = sum(abs2, trial_f)

        # Step quality = residual change / predicted residual change
        rho = (T_resid - C_resid) / (P_resid - C_resid)
        if rho > min_step_quality
            # Update x, 𝐹, and residual
            x .= xnew
            value!(df, x)
            C_resid = T_resid

            # Increase trust region radius
            if rho > good_step_quality
                λ = max(λᵣ * λ, λₘ)
            end
        else
            # Decrease trust region radius
            λ = min(λᵢ * λ, Λₘ)
        end

        # Increase the iteration
        iter += 1

        # Check convergence criteria:
        # 1. Small gradient: norm(𝐽ᵀ * 𝐹, Inf) < g_tol
        if norm(𝐽' * 𝐹, Inf) < g_tol
            g_converged = true
        end
        # 2. Small step size: norm(δx) < x_tol
        if norm(δx) < x_tol * (x_tol + norm(x))
            x_converged = true
        end
        converged = g_converged | x_converged
    end

    # Return the results
    OptimizationResults(
        x₀,          # x₀
        x,           # minimizer
        C_resid,     # minimum (residual)
        iter,        # iterations
        !converged,  # iteration_converged
        x_converged, # x_converged
        g_converged, # g_converged
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
