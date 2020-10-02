# printing for static taylor
function coeffstring(t::STaylor1, i, variable=:t)
    if i == 1  # order 0
        return string(t.coeffs[i])
    end
    if i == 2  # order 1
        return string(t.coeffs[i], variable)
    end
    return string(t.coeffs[i], variable, "^", i-1)
end

function print_taylor(io::IO, t::STaylor1, variable=:t)
    print(io, "(" * join([coeffstring(t, i, variable) for i in 1:length(t.coeffs)], " + ") * ")")
end

Base.show(io::IO, t::STaylor1) = print_taylor(io, t)
function Base.show(io::IO, t::STaylor1{N,T}) where {N, T<:STaylor1}
    print_taylor(io, t, :t)
end
