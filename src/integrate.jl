module Integrate
export romberg, int2d

function romberg(f::Function, a::Real, b::Real, e::Real)::Real
    table = Dict()
    res = table[1,1] = (b - a) * (f(a) + f(b)) / 2
    for j in 2 : 50
        hⱼ = (b - a) / 2 ^ (j - 1)
        table[j,1] = 0.5 * table[j-1,1] + hⱼ * sum(f(a + (2i - 1)hⱼ) for i in 1 : 2^(j - 2))
        for k in 2 : j
            table[j,k] = (4^(k - 1) * table[j,k-1] - table[j-1,k-1]) / (4^(k-1) - 1)
        end

        res = table[j,j]
        if abs(table[j,j] - table[j-1,j-1]) < e
            break
        end
    end

    res
end

"""
A naive two dimensional integral method.
"""
function int2d(f::Function, (a, b, m)::Tuple{<:Real, <:Real, Int}, (c, d, n)::Tuple{<:Real, <:Real, Int})::Real
    h = (b - a) / m
    k = (d - c) / n
    res = 0.0
    for i in 0 : m
        for j in 0 : n
            res += f(a + i * h, c + j * k) * if i in [0, m] && j in [0, n]
                0.25
            elseif i in [0, m] || j in [0, n]
                0.5
            else
                1.0
            end
        end
    end
    h * k * res
end

end