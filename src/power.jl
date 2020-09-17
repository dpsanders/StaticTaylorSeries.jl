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

@generated function inverse(a::STaylor1{N,T}) where {N,T<:Real}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:(2*N)])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    exout = :(($(syms[1]),))
    for i = 2:N
        push!(exout.args, syms[i])
    end

    count = 1
    for n = 1:(N - 1)
        ex_calc.args[2*count - 1] = :(syms[n + 1] = zdivfpown[n - 1]/n)
        ex_calc.args[2*count] = :(zdivfpown *= zdivf)
        count += 1
    end

    return quote
               Base.@_inline_meta
               if a[0] != zero(T)
                   throw(ArgumentError(
                   """
                   Evaluation of Taylor1 series at 0 is non-zero. For high accuracy, revert
                   a Taylor1 series with first coefficient 0 and re-expand about f(0).
                   """))
               end
               z = copy(a)
               zdivf = z/a
               zdivfpown = zdivf
               S = eltype(zdivf)
               $ex_calc
               return STaylor1{N,T}($exout)
            end
end


@generated function sqrt(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    start_expr = quote
                    Base.@_inline_meta
                    l0nz = findfirst(a)
                    aux = zero(T)
                    if l0nz < 0
                        return
                    elseif l0nz%2 == 1
                        throw(ArgumentError(
                        """First non-vanishing Taylor1 coefficient must correspond
                        to an **even power** in order to expand `sqrt` around 0."""))
                    end
                    lnull = div(l0nz, 2)
                 end

    for k = (lnull + 1):(N-1)

        pre_loop = quote
            if k == lnull
                @inbounds $(syms[k]) = sqrt(aa[2*lnull])
                return continue
            end
            kodd = (k - lnull)%2
            kend = div(k - lnull - 2 + kodd, 2)
            imax = min(lnull + kend, N - 1)
            imin = max(lnull + 1, k + lnull - N + 1)
            imin ≤ imax && (@inbounds $(syms[k]) = $(syms[imin])*$(syms[k + lnull - imin]))
        end

        loop = quote
        end

        post_loop = quote
        end


        kT = convert(T,k)
        sym = syms[k+1]
        ex_line = :($kT * a[$k] * $(syms[1]))
        @inbounds for i = 1:k-1
            ex_line = :($ex_line + $(kT-i) * a[$(k-i)] * $(syms[i+1]))
        end
        ex_line = :(($ex_line)/$kT)
        ex_line = :($sym = $ex_line)
        ex_calc.args[k+1] = ex_line
    end

    exout = :(($(syms[1]),))
    for i = 2:N
        push!(exout.args, syms[i])
    end

    return quote
        $start_expr
        $compute_expr
        return STaylor1{N,T}($exout)
    end
end
