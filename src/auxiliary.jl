# getindex STaylor1
getindex(a::STaylor1, n::Int) = a.coeffs[n+1]
getindex(a::STaylor1, u::UnitRange{Int}) = a.coeffs[u .+ 1]
getindex(a::STaylor1, c::Colon) = a.coeffs[c]
getindex(a::STaylor1, u::StepRange{Int,Int}) = a.coeffs[u .+ 1]

@inline iterate(a::STaylor1{N,T}, state=0) where {N, T<:Number} = state > N-1 ? nothing : (a.coeffs[state+1], state+1)
@inline firstindex(a::STaylor1) = 0
@inline lastindex(a::STaylor1{N,T}) where {N, T<:Number} = N-1
@inline eachindex(s::STaylor1{N,T}) where {N, T<:Number} = UnitRange(0, N-1)
@inline size(s::STaylor1{N,T}) where {N, T<:Number} = N
@inline length(s::STaylor1{N,T}) where {N, T<:Number} = N
@inline get_order(s::STaylor1{N,T}) where {N, T<:Number} = N - 1
@inline eltype(s::STaylor1{N,T}) where {N, T<:Number} = T
@inline axes(a::STaylor1) = ()

function Base.findfirst(a::STaylor1{N,T}) where {N, T<:Number}
    first = findfirst(x->!iszero(x), a.coeffs)
    isa(first, Nothing) && return -1
    return first-1
end
# Finds the last non-zero entry; extended to Taylor1
function Base.findlast(a::STaylor1{N,T}) where {N, T<:Number}
    last = findlast(x->!iszero(x), a.coeffs)
    isa(last, Nothing) && return -1
    return last-1
end
constant_term(a::STaylor1) = a[0]
