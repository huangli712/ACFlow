
function read_data!(ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    niw = 64
    #
    open("giw.data", "r") do fin
        for i = 1:niw
            arr = parse.(F64, line_to_array(fin))
            push!(ω.grid, arr[1])
            push!(𝐺.value, arr[2] + arr[3] * im)
        end
    end
    #
    open("err.data", "r") do fin
        for i = 1:niw
            arr = parse.(F64, line_to_array(fin))
            @assert ω.grid[i] == arr[1]
            push!(𝐺.error, arr[2] + arr[3] * im)
        end
    end
end

function read_data!(τ::ImaginaryTimeGrid, 𝐺::GreenData)

end