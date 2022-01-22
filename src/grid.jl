
function Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.ngrid
    return fg.ω[ind]
end

function make_grid(fg::FermionicMatsubaraGrid)
    for n in eachindex(fg.ω)
        fg.ω[n] = (2 * n - 1) * π / fg.β
    end
end