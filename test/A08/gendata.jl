#!/usr/bin/env julia

using Random
using Printf

# Numerical integration by the composite trapezoidal rule
function trapz(x, y, linear::Bool = false)
    # For linear mesh
    if linear
        h = x[2] - x[1]
        value = y[1] + y[end] + 2.0 * sum(y[2:end-1])
        value = h * value / 2.0
    # For non-equidistant mesh
    else
        len = length(x)
        value = 0.0
        for i = 1:len-1
            value = value + (y[i] + y[i+1]) * (x[i+1] - x[i])
        end
        value = value / 2.0
    end

    return value
end

# Setup parameters
wmin = -5.0  # Left boundary
wmax = +5.0  # Right boundary
nmesh = 5001 # Number of real-frequency points
niw  = 20    # Number of Matsubara frequencies
beta = 40.0  # Inverse temperature

# Real frequency mesh
w_real = collect(LinRange(wmin, wmax, nmesh))

# Initial spectral function
spec_real1 = similar(w_real)
@. spec_real1  = 0.5 * exp(-(w_real - 1.0) ^ 2.0 / (2.0 * 0.2 ^ 2.0)) / (0.2 * sqrt(2.0 * π))
@. spec_real1 += 0.5 * exp(-(w_real - 2.0) ^ 2.0 / (2.0 * 0.7 ^ 2.0)) / (0.7 * sqrt(2.0 * π))
#
spec_real2 = similar(w_real)
@. spec_real2  = 0.5 * exp(-(w_real + 1.0) ^ 2.0 / (2.0 * 0.25^ 2.0)) / (0.25* sqrt(2.0 * π))
@. spec_real2 += 0.5 * exp(-(w_real + 2.1) ^ 2.0 / (2.0 * 0.6 ^ 2.0)) / (0.6 * sqrt(2.0 * π))
#
spec_matrix = zeros(Float64, (2,2,nmesh))
spec_matrix[1,1,:] .= spec_real1
spec_matrix[2,2,:] .= spec_real2

# Rotate spectral function to generate non-diagonal element
#
# Rotation angle
rot_ang = 0.1
#
# Rotation matrix
rot_mat = [cos(rot_ang) sin(rot_ang); -sin(rot_ang) cos(rot_ang)]
T_rot_mat = rot_mat'
#
# Get final spectral function
true_spec = zeros(Float64, (2,2,nmesh))
for i = 1:2
    for l = 1:2
        for w = 1:nmesh
            for j = 1:2
                for k = 1:2
                    true_spec[i,l,w] = true_spec[i,l,w] +
                        rot_mat[i,j] * spec_matrix[j,k,w] * T_rot_mat[k,l]
                end
            end
        end
    end
end

# Matsubara frequency mesh
iw = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Kernel function
kernel = 1.0 ./ (im * reshape(iw, (niw,1)) .- reshape(w_real, (1,nmesh)))

# Build green's function
KA = reshape(kernel, (1,1,niw,nmesh)) .* reshape(true_spec, (2,2,1,nmesh))
giw = zeros(ComplexF64, (2,2,niw))
for i = 1:2
    for j = 1:2
        for w = 1:niw
            giw[i,j,w] = trapz(w_real, KA[i,j,w,:])
        end
    end
end

# Build error
err = 1e-5

# Write green's function
open("green.11.data", "w") do fout
    for i = 1:niw
        z = giw[1,1,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
open("green.12.data", "w") do fout
    for i = 1:niw
        z = giw[1,2,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
open("green.22.data", "w") do fout
    for i = 1:niw
        z = giw[2,2,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end

# Write spectral function
open("exact.11.data", "w") do fout
    for i in eachindex(w_real)
        @printf(fout, "%20.16f %20.16f\n", w_real[i], true_spec[1,1,i])
    end
end
open("exact.12.data", "w") do fout
    for i in eachindex(w_real)
        @printf(fout, "%20.16f %20.16f\n", w_real[i], true_spec[1,2,i])
    end
end
open("exact.22.data", "w") do fout
    for i in eachindex(w_real)
        @printf(fout, "%20.16f %20.16f\n", w_real[i], true_spec[2,2,i])
    end
end
