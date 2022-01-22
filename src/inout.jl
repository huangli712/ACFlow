
function read_data()
    finput = get_c("finput")
    ngrid = get_c("ngrid")
    grid = get_c("grid")
    @show finput, ngrid, grid

    if grid == "matsubara"
        return read_matsubara_data(finput, ngrid)
    else
        return read_imaginary_time_data(finput, ngrid)
    end
end

function read_imaginary_time_data(finput::String, ngrid::I64)
end

function read_matsubara_data(finput::String, ngrid::I64)
    mesh = zeros(F64, ngrid)
    value = zeros(C64, ngrid)
    error = zeros(C64, ngrid)

    open(finput, "r") do fin
        for i = 1:ngrid
            arr = parse.(F64, line_to_array(fin)[1:3])
            mesh[i] = arr[1]
            value[i] = arr[2] + im * arr[3]
        end
    end

    return RawData(mesh, value, error)
end