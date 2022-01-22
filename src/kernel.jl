function make_kernel(um::UniformMesh, fg::FermionicMatsubaraGrid)
    niw = fg.ngrid
    nw = um.nmesh

    kernel = zeros(C64, niw, nw)
    for i = 1:nw
        for j = 1:niw
            kernel[j,i] = 1.0 / (im * fg[j] - um[i])
        end
    end

    return kernel
end