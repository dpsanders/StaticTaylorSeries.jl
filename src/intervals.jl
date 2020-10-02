using .IntervalArithmetic

function evaluate(a::STaylor1{N,T}, dx::Interval) where {N, T <: Number}
    order = N - 1
    uno = one(dx)
    dx2 = dx^2
    if iseven(order)
        kend = order-2
        @inbounds sum_even = a[end]*uno
        @inbounds sum_odd = a[end-1]*zero(dx)
    else
        kend = order-3
        @inbounds sum_odd = a[end]*uno
        @inbounds sum_even = a[end-1]*uno
    end
    @inbounds for k = kend:-2:0
        sum_odd = sum_odd*dx2 + a[k + 1]
        sum_even = sum_even*dx2 + a[k]
    end
    return sum_even + sum_odd*dx
end

normalize_taylor(a::STaylor1{N,T}, I::Interval{T}, symI::Bool=true) where {N, T <: Number} =
    _normalize(a, I, Val(symI))

function _normalize(a::STaylor1{N,T}, I::Interval{T}, ::Val{true}) where {N, T <: Number}
    t = STaylor1(zero(T), Val{N}())
    tnew = mid(I) + t*radius(I)
    return a(tnew)
end

function _normalize(a::STaylor1{N,T}, I::Interval{T}, ::Val{false}) where {N, T <: Number}
    t = STaylor1(zero(promote_type(S,T)), Val{N}())
    tnew = inf(I) + t*diam(I)
    return a(tnew)
end
