function make_model()
    model = get_c["model"]
    if model == "flat"
        return make_flat_model()
    else
        return make_gaussian_model()
    end
end

function make_flat_model()
end

function make_gaussian_model()
end