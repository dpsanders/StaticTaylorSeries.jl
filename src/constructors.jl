######################### STaylor1
"""
    STaylor1{N,T<:Number} <: AbstractSeries{T}

DataType for polynomial expansions in one independent variable.

**Fields:**

- `coeffs :: NTuple{N,T}` Expansion coefficients; the ``i``-th
    component is the coefficient of degree ``i-1`` of the expansion.

Note that `STaylor1` variables are callable. For more information, see
[`evaluate`](@ref).
"""
struct STaylor1{N,T<:Number} <: AbstractSeries{T}
    coeffs::NTuple{N,T}
    function STaylor1{N,T}(coeffs::NTuple{N,T}) where {N, T <: Number}
        new(coeffs)
    end
end
function STaylor1(coeffs::NTuple{N,T}) where {N, T <: Number}
    STaylor1{N,T}(coeffs)
end

## Outer constructors ##

"""
    STaylor1(x::T, v::Val{N})

Shortcut to define the independent variable of a `STaylor1{N,T}` polynomial of
given `N` with constant term equal to `x`.
"""
@generated function STaylor1(x::T, v::Val{N}) where {N,T<:Number}
    y = Any[:(zero($T)) for i=1:N]
    tup = :((x,))
    push!(tup.args, y...)
    return quote
        Base.@_inline_meta
        STaylor1{(N+1),T}($tup)
    end
end
function STaylor1(coeffs::Vector{T}, l::Val{L}, v::Val{N}) where {L,N,T<:Number}
    STaylor1{(N+1),T}(ntuple(i -> (i > L+1) ? coeffs[i] : zero(T),  N+1))
end
@inline function STaylor1(coeffs::Vector{T}) where {T<:Number}
    STaylor1{length(coeffs),T}(tuple(coeffs...))
end
@inline STaylor1(x::STaylor1{N,T}) where {N,T<:Number} = x
