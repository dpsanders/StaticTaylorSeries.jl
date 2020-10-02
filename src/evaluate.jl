function evaluate(a::STaylor1{N,T}, dx::T) where {N, T<:Number}
    @inbounds suma = a[N-1]
    @inbounds for k in (N-1):-1:0
        suma = suma*dx + a[k]
    end
    suma
end
#=
function evaluate(a::STaylor1{N,T}, dx::S) where {N, T<:Number, S<:Number}
    suma = a[N-1]*one(dx)
    @inbounds for k in (N-1):-1:0
        suma = suma*dx + a[k]
    end
    suma
end
=#
evaluate(a::STaylor1{N,T}) where {N, T<:Number} = a[0]

evaluate(x::Union{Array{STaylor1{N,T}}, SubArray{STaylor1{N,T}}}, δt::S) where
        {N, T<:Number, S<:Number} = evaluate.(x, δt)
evaluate(a::Union{Array{STaylor1{N,T}}, SubArray{STaylor1{N,T}}}) where
        {N, T<:Number} = evaluate.(a, zero(T))

(p::STaylor1)(x) = evaluate(p, x)
(p::STaylor1)() = evaluate(p)

(p::Array{STaylor1{N,T}})(x) where {N,T<:Number} = evaluate.(p, x)
(p::SubArray{STaylor1{N,T}})(x) where {N,T<:Number} = evaluate.(p, x)
(p::Array{STaylor1{N,T}})() where {N,T<:Number} = evaluate.(p)
(p::SubArray{STaylor1{N,T}})() where {N,T<:Number} = evaluate.(p)
