function ^(a::STaylor1{N,T}, n::Integer) where {N,T<:Real}
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return TaylorSeries.square(a)
    #n < 0 && return a^float(n)
    return power_by_squaring(a, n)
end

^(a::STaylor1{N,T}, b::STaylor1{N,T}) where {N,T<:Number} = exp(b*log(a))

function power_by_squaring(x::STaylor1{N,T}, p::Integer) where {N,T<:Number}
    p == 1 && return x
    p == 0 && return one(x)
    p == 2 && return square(x)
    t = trailing_zeros(p) + 1
    p >>= t

    while (t -= 1) > 0
        x = square(x)
    end

    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) ≥ 0
            x = square(x)
        end
        y *= x
    end

    return y
end


@generated function square(a::STaylor1{N,T}) where {N, T<:Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(syms[1]) = a[0]^2)
    ex_calc.args[1] = ex_line

    for k in 1:(N-1)
        kodd = k%2
        kend = div(k - 2 + kodd, 2)
        ex_line = :(a[0] * a[$k])
        @inbounds for i = 1:kend
            ex_line = :($ex_line + a[$i] * a[$(k-i)])
        end
        ex_line = :(2.0*($ex_line)) # float(2)* TODO: ADD BACK IN
        if kodd !== 1
            ex_line = :($ex_line +a[$(div(k,2))]^2)
        end
        ex_line = :($(syms[k+1]) = $ex_line)
        ex_calc.args[k+1] = ex_line
    end

    exout = :(($(syms[1]),))
    for i = 2:N
        push!(exout.args, syms[i])
    end
    return quote
               Base.@_inline_meta
               $ex_calc
               return STaylor1{N,T}($exout)
            end
end
