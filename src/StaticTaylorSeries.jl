module StaticTaylorSeries

export StaticTaylor

"""
Static Taylor series with variable V, length N (i.e. order N-1) and element type T.
"""
struct StaticTaylor{N,T}
    coeffs::NTuple{N,T}
end


# function print_taylor(io::IO, t::StaticTaylor, variable=:x)



function coeffstring(t::StaticTaylor, i, variable=:x)

    if i == 1  # order 0
        return string(t.coeffs[i])

    elseif i == 2  # order 1
        return string(t.coeffs[i], variable)

    else
        return string(t.coeffs[i], variable, "^", i-1)
    end
end


function print_taylor(io::IO, t::StaticTaylor, variable=:x)

    print(io, "(" * join([coeffstring(t, i, variable) for i in 1:length(t.coeffs)], " + ") * ")")

end

function Base.show(io::IO, t::StaticTaylor)
    print_taylor(io, t)
end

function Base.show(io::IO, t::StaticTaylor{N,T}) where {N, T<:StaticTaylor}
    print_taylor(io, t, :y)
end

#StaticTaylor(iterable...) = StaticTaylor(SVector(iterable...))

StaticTaylor{N}(v::NTuple{N,T}) where {N,T} = StaticTaylor(v)

StaticTaylor{N}(iterable) where {N} = StaticTaylor{N}(SVector{N}(iterable))

StaticTaylor{N}(iterable...) where {N} = StaticTaylor{N}(iterable)

# StaticTaylor(v) = StaticTaylor(v)

StaticTaylor(iterable...) = StaticTaylor{length(iterable)}(iterable)

import Base:getindex, length, eltype

getindex(s::StaticTaylor, i::Integer) = s.coeffs[i+1]

length(s::StaticTaylor{N,T}) where {N,T} = N

eltype(s::StaticTaylor{N,T}) where {N,T} = T


import Base: +, -, * ,^

function +(s::StaticTaylor{N,T}, t::StaticTaylor{N,T}) where {N,T}
    return StaticTaylor(s.coeffs .+ t.coeffs)
end

function +(s::StaticTaylor{N,T}, α::Real) where {N,T}
    StaticTaylor{N,T}(ntuple(i -> i == 1 ? s.coeffs[1] + α : s.coeffs[i], N))
end

+(α::Real, s::StaticTaylor) == s + α

-(s::StaticTaylor{N,T}) where {N,T} = StaticTaylor(.-(s.coeffs))

function -(s::StaticTaylor{N,T}, t::StaticTaylor{N,T}) where {N,T}
    return StaticTaylor(s.coeffs .- t.coeffs)
end

-(s::StaticTaylor{N,T}, α::Real) = s + (-α)

-(α::Real, s::StaticTaylor{N,T}) = -(s - a)





Base.literal_pow(::typeof(^), x::StaticTaylor, ::Val{p}) where {p} = x^p

^(x::StaticTaylor, n::Integer) = Base.power_by_squaring(x, n)



# function *(s::StaticTaylor{N,T}, t::StaticTaylor{N,T}) where {N,T}
#     v = SVector(ntuple(k->sum(s[i]*t[k-1-i] for i in 0:k-1), Val(N)))
#     return StaticTaylor(v)
# end

# The following is modified from StaticUnivariatePolynomials.jl
@generated function Base.:*(p1::StaticTaylor{N,T}, p2::StaticTaylor{N,T}) where {N, T}
    exprs = Any[nothing for i in 1:N]
    for i in 0 : N-1   # order is N-1
        for j in 0 : N-1
            k = i + j + 1  # setindex does not have offset

            if k > N
                continue
            end

            if exprs[k] === nothing
                exprs[k] = :(p1[$i] * p2[$j])
            else
                exprs[k] = :(muladd(p1[$i], p2[$j], $(exprs[k])))
            end
        end
    end

    # Core.println("Generated code with N=$N:")
    # Core.println(exprs)
    # Core.println()

    return quote
        Base.@_inline_meta
        StaticTaylor{N,T}(tuple($(exprs...)))
    end
end

*(s::StaticTaylor, α::Real) = StaticTaylor(α .* s.coeffs)
*(α::Real, s::StaticTaylor) = s * α

/(s::StaticTaylor, α::Real) = StaticTaylor(s.coeffs ./ α)

end # module
