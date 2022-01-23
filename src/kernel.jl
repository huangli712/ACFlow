function make_kernel(m::AbstractMesh, fg::FermionicMatsubaraGrid)
    niw = fg.nfreq
    nw = m.nmesh

    kernel = zeros(C64, niw, nw)
    for i = 1:nw
        for j = 1:niw
            kernel[j,i] = 1.0 / (im * fg[j] - m[i])
        end
    end

    return kernel
end