

function make_kernel(am::AbstractMesh, fg::FermionicImaginaryTimeGrid)
end

function make_kernel(am::AbstractMesh, fg::FermionicMatsubaraGrid)
    nfreq = fg.nfreq
    nmesh = am.nmesh

    kernel = zeros(C64, nfreq, nmesh)
    for i = 1:nmesh
        for j = 1:nfreq
            kernel[j,i] = 1.0 / (im * fg[j] - am[i])
        end
    end

    return kernel
end

function make_kernel(am::AbstractMesh, bg::BosonicImaginaryTimeGrid)
end

function make_kernel(am::AbstractMesh, bg::BosonicMatsubaraGrid)
end