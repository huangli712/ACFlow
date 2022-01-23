
function make_model(m::AbstractMesh)
    model = get_c("model")
    if model == "flat"
        return make_flat_model(m)
    else
        return make_gaussian_model(m)
    end
end

function make_flat_model(um::UniformMesh)
    model = ones(F64, um.nmesh)
    norm = trapz(um.mesh, model)
    model = model ./ norm
    return model
end

function make_flat_model(num::NonUniformMesh)
end

function make_gaussian_model(um::UniformMesh)
end

function make_gaussian_model(num::NonUniformMesh)
end