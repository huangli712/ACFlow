function read_time_data(finput::String, ngrid::I64)
end

function read_freq_data(finput::String, ngrid::I64)
    mesh = zeros(F64, ngrid)
    value = zeros(C64, ngrid)
    error = zeros(C64, ngrid)

    open(finput, "r") do fin
        for i = 1:ngrid
            arr = parse.(F64, line_to_array(fin)[1:3])
            mesh[i] = arr[1]
            value[i] = arr[2] + im * arr[3]
            error[i] = 0.0001 + 0.0001im
        end
    end

    return RawData(mesh, value, error)
end