# Conversion for STaylor1
function convert(::Type{STaylor1{N,Rational{T}}}, a::STaylor1{N,S}) where {N,T<:Integer, S<:AbstractFloat}
    STaylor1{N,T}(rationalize.(a[:]))
end
function convert(::Type{STaylor1{N,T}}, b::Array{T,1}) where {N,T<:Number}
    @assert N == length(b)
    STaylor1{N,T}(b)
end
function convert(::Type{STaylor1{N,T}}, b::Array{S,1}) where {N,T<:Number, S<:Number}
    @assert N == length(b)
    STaylor1{N,T}(convert(Array{T,1},b))
end
convert(::Type{STaylor1{N,T}}, a::STaylor1{N,T}) where {N,T<:Number} = a
convert(::Type{STaylor1{N,T}}, b::S)  where {N, T<:Number, S<:Number} = STaylor1(convert(T,b), Val(N))
convert(::Type{STaylor1{N,T}}, b::T)  where {N, T<:Number} = STaylor1(b, Val(N))

promote_rule(::Type{STaylor1{N,T}}, ::Type{STaylor1{N,T}}) where {N, T<:Number} = STaylor1{N,T}
promote_rule(::Type{STaylor1{N,T}}, ::Type{STaylor1{N,S}}) where {N, T<:Number, S<:Number} = STaylor1{N, promote_type(T,S)}
promote_rule(::Type{STaylor1{N,T}}, ::Type{Array{T,1}}) where {N, T<:Number} = STaylor1{N,T}
promote_rule(::Type{STaylor1{N,T}}, ::Type{Array{S,1}}) where {N, T<:Number, S<:Number} = STaylor1{N,promote_type(T,S)}
promote_rule(::Type{STaylor1{N,T}}, ::Type{T}) where {N, T<:Number} = STaylor1{N,T}
promote_rule(::Type{STaylor1{N,T}}, ::Type{S}) where {N, T<:Number, S<:Number} = STaylor1{N,promote_type(T,S)}
