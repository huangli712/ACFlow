function Base.getindex(um::UniformMesh, ind::I64)
    @assert 1 ≤ ind ≤ um.nmesh
    return um.mesh[ind]
end

function make_uniform_mesh()
    nmesh = get_c["nmesh"]
    wmax = get_c["wmax"]
    wmin = get_c["wmin"]

    mesh = collect(LinRange(wmin, wmax, nmesh))
    weight = (mesh[2:end] + mesh[1:end-1]) / 2.0
    pushfirst!(weight, mesh[1])
    push!(weight, mesh[end])
    weight = diff(weight)

    return UniformMesh(nmesh, wmax, wmin, mesh, weight)
end

function make_non_uniform_mesh()
end