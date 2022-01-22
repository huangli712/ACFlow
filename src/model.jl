function make_model(um::UniformMesh)
    model = get_c("model")
    if model == "flat"
        return make_flat_model(um)
    else
        return make_gaussian_model(um)
    end
end

function make_flat_model(um::UniformMesh)
    model = ones(F64, um.nmesh)
    norm = trapz(um.mesh, model)
    model = model ./ norm

    return model
end

function make_gaussian_model(um::UniformMesh)
end