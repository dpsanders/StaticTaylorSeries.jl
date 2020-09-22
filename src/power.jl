function ^(a::STaylor1{N,T}, n::Integer) where {N,T<:Real}
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return square(a)
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

@generated function ^(a::STaylor1{N,T}, r::S) where {N, T<:Number, S<:Real}

    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]
    ctuple = Expr(:tuple)
    for i = 1:N
        push!(ctuple.args, syms[i])
    end

    for i = 1:N
        push!(ex_calc.args, :($(syms[i]) = zero(T)))
    end

    expr_quote = quote
        iszero(r) && return one(a)
        r == 1 && return a
        r == 2 && return square(a)
        r == 1/2 && return sqrt(a)
        $ex_calc
    end

    c = STaylor1(zero(T), Val{N}())
    for k = 0:(N - 1)
        symk = syms[k + 1]
        temp_quote = quote
            # First non-zero coefficient
            l0 = findfirst(a)
            if l0 < 0
                $symk = zero(T)
            else
                # The first non-zero coefficient of the result; must be integer
                !isinteger(r*l0) && throw(ArgumentError(
                    """The 0th order Taylor1 coefficient must be non-zero
                    to raise the Taylor1 polynomial to a non-integer exponent."""))
                lnull = trunc(Int, r*l0)
                kprime = k - lnull
                if (kprime < 0) || (lnull > a.order)
                    $symk = zero(T)
                else
                    # Relevant for positive integer r, to avoid round-off errors
                    if isinteger(r) && (k > r*findlast(a))
                        $symk = zero(T)
                    else
                        if k == lnull
                            $symk = a[l0]^r
                        else
                            # The recursion formula
                            if l0 + kprime ≤ (N - 1)
                                tup_in = $ctuple
                                $symk = r*kprime*tup_sel(lnull, tup_in)*a[l0 + kprime]
                            else
                                $symk = zero(T)
                            end
                            for i = 1:(k - lnull - 1)
                                if !((i + lnull) > (N - 1) || (l0 + kprime - i > (N - 1)))
                                    aux = r*(kprime - i) - i
                                    tup_in = $ctuple
                                    $symk += aux*tup_sel(i + lnull, tup_in)*a[l0 + kprime - i]
                                end
                            end
                            $symk /= kprime*a[l0]
                        end
                    end
                end
            end
        end
        expr_quote = quote
            $expr_quote
            if r == 0
                # one(c)
            elseif r == 1
                # DO NOTHING
            elseif r == 2
                # square(c)
            elseif r == 0.5
                # sqrt(c)
            else
                $temp_quote
            end
        end
    end

    exout = :(($(syms[1]),))
    for i = 2:N
        push!(exout.args, syms[i])
    end

    return quote
               Base.@_inline_meta
               $expr_quote
               return STaylor1{N,T}($exout)
            end
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

function tup_sel(i, vargs)
    return vargs[i+1]
end

@generated function sqrt(a::STaylor1{N,T}) where {N,T<:Number}

    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]
    ctuple = Expr(:tuple)
    for i = 1:N
        push!(ctuple.args, syms[i])
    end

    # First non-zero coefficient
    expr_quote = quote
        l0nz = findfirst(a)
        aux = zero(T)
        if l0nz < 0
            return zero(STaylor1{N,T})
        elseif l0nz%2 == 1 # l0nz must be pair
            throw(ArgumentError(
            """First non-vanishing Taylor1 coefficient must correspond
            to an **even power** in order to expand `sqrt` around 0."""))
        end

        # The last l0nz coefficients are set to zero.
        lnull = div(l0nz, 2)
    end

    for i = 1:N
        push!(ex_calc.args, :($(syms[i]) = zero(T)))
    end


    for i = 1:N
        switch_expr = :((lnull == $(i-1)) && ($(syms[i]) = sqrt(a[l0nz])))
        expr_quote = quote
            $expr_quote
            $switch_expr
        end
    end

    for k = 0:(N - 1)
        symk = syms[k + 1]
        temp_expr = quote
            if $k >= lnull + 1
                if $k == lnull
                    $symk = sqrt(a[2*lnull])
                else
                    kodd = ($k - lnull)%2
                    kend = div($k - lnull - 2 + kodd, 2)
                    imax = min(lnull + kend, N - 1)
                    imin = max(lnull + 1, $k + lnull - (N - 1))
                    if imin ≤ imax
                        tup_in = $ctuple
                        $symk = tup_sel(imin, tup_in)*tup_sel($k + lnull - imin, tup_in)
                    end
                    for i = (imin + 1):imax
                        tup_in = $ctuple
                        $symk += tup_sel(i, tup_in)*tup_sel($k + lnull - i, tup_in)
                    end
                    if $k + lnull ≤ (N - 1)
                        aux = a[$k + lnull] - 2*$symk
                    else
                        aux = -2*$symk
                    end
                    tup_in = $ctuple
                    if kodd == 0
                        aux -= tup_sel(kend + lnull + 1, tup_in)^2
                    end
                    $symk = aux/(2*tup_sel(lnull, tup_in))
                end
            end
        end
        expr_quote = quote
            $expr_quote
            $temp_expr
        end
    end

    exout = :(($(syms[1]),))
    for i = 2:N
        push!(exout.args, syms[i])
    end
    return quote
               Base.@_inline_meta
               $ex_calc
               $expr_quote
               return STaylor1{N,T}($exout)
            end
end
