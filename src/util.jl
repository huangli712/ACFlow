@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

#function myfun(a, b)
#    #return a ^ 2.0 + b * a + 1.0
#    return a ^ 3.0 - b
#end

#s = secant(myfun, 1.0, 8.0)
#println(s)

function secant(func, x0, args)
    eps = 1.0e-4
    maxiter = 50
    tol = 1.48e-8
    funcalls = 0
    p0 = 1.0 * x0
    p1 = x0 * (1.0 + eps)
    if p1 ≥ 0.0
        p1 = p1 + eps
    else
        p1 = p1 - eps
    end

    q0 = func(p0, args)
    funcalls = funcalls + 1
    q1 = func(p1, args)
    funcalls = funcalls + 1

    if abs(q1) < abs(q0)
        p0, p1 = p1, p0
        q0, q1 = q1, q0
    end

    for itr = 1:maxiter
        if q1 == q0
            if p1 != p0
                error("tolerance is reached!")
            end
            p = (p1 + p0) / 2.0
            return p
        else
            if abs(q1) > abs(q0)
                p = (-q0 / q1 * p1 + p0) / (1 - q0 / q1)
            else
                p = (-q1 / q0 * p0 + p1) / (1 - q1 / q0)
            end
        end

        if abs(p - p1) < tol
            return p
        end

        p0, q0 = p1, q1
        p1 = p
        q1 = func(p1, args)
        funcalls = funcalls + 1
    end
end

function trapz(x, y)
    #@show x
    #@show y
    h = x[2] - x[1]
    sum = 0.0
    for i = 2:length(x)-1
        sum = sum + y[i]
    end
    value = (h / 2.0) * (y[1] + y[end] + 2.0 * sum)
        
    return value
end

function new_trapz(x, y)
    len = length(x)
    value = 0.0
    for i = 1:len-1
        value = value + (y[i] + y[i+1]) * (x[i+1] - x[i]) / 2.0
    end
    return value
end