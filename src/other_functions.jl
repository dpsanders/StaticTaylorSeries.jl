for f in (:real, :imag, :conj)
    @eval ($f)(a::STaylor1{N,T}) where {N,T<:Number} = STaylor1{N,T}(($f).(a.coeffs))
end

adjoint(a::STaylor1) = conj(a)
isinf(a::STaylor1) = any(isinf.(a.coeffs))
isnan(a::STaylor1) = any(isnan.(a.coeffs))

function abs(a::STaylor1{N,T}) where {N,T<:Real}
    if a[0] > zero(T)
        return a
    elseif a[0] < zero(T)
        return -a
    else
        throw(ArgumentError(
        """The 0th order Taylor1 coefficient must be non-zero
        (abs(x) is not differentiable at x=0)."""))
    end
end

abs2(a::STaylor1{N,T}) where {N,T<:Real} = a^2

deg2rad(z::STaylor1{N, T}) where {N, T<:AbstractFloat} = z * (convert(T, pi) / 180)
deg2rad(z::STaylor1{N, T}) where {N, T<:Real} = z * (convert(float(T), pi) / 180)

rad2deg(z::STaylor1{N, T}) where {N, T<:AbstractFloat} = z * (180 / convert(T, pi))
rad2deg(z::STaylor1{N, T}) where {N, T<:Real} = z * (180 / convert(float(T), pi))
