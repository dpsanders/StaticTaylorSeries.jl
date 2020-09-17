# Functions for STaylor1
@generated function exp(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(syms[1]) = exp(a[0]))
    ex_calc.args[1] = ex_line

    for k in 1:(N-1)
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
               Base.@_inline_meta
               $ex_calc
               return STaylor1{N,T}($exout)
            end
end

@generated function log(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    (N >= 1) && (ex_calc.args[1] = :($(syms[1]) = log(constant_term(a))))
    (N >= 2) && (ex_calc.args[2] = :($(syms[2]) = a[1]/constant_term(a)))

    for k in 2:(N-1)
        ex_line = :($(k-1)*a[1]*$(syms[k]))
        @inbounds for i = 2:k-1
            ex_line = :($ex_line + $(k-i)*a[$i] * $(syms[k+1-i]))
        end
        ex_line = :((a[$k] - ($ex_line)/$(convert(T,k)))/constant_term(a))
        ex_line = :($(syms[k+1]) = $ex_line)
        ex_calc.args[k+1] = ex_line
    end

    exout = :(($(syms[1]),))
    for i = 2:N
        push!(exout.args, syms[i])
    end

    return quote
               Base.@_inline_meta
               iszero(constant_term(a)) && throw(ArgumentError("""
                       The 0-th order `STaylor1` coefficient must be non-zero
                       in order to expand `log` around 0.
                       """))
               $ex_calc
               return STaylor1{N,T}($exout)
            end
end

sin(a::STaylor1{N,T}) where {N, T <: Number} = sincos(a)[1]
cos(a::STaylor1{N,T}) where {N, T <: Number} = sincos(a)[2]
@generated function sincos(a::STaylor1{N,T}) where {N, T <: Number}

    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:(2*N)])

    syms_s = Symbol[Symbol("c$i") for i in 1:N]
    syms_c = Symbol[Symbol("c2$i") for i in 1:N]

    ex_line_s = :($(syms_s[1]) = sin(a[0]))
    ex_line_c = :($(syms_c[1]) = cos(a[0]))
    ex_calc.args[1] = ex_line_s
    ex_calc.args[2] = ex_line_c

    for k = 1:(N - 1)
        ex_line_s = :(a[1]*$(syms_c[k]))
        ex_line_c = :(-a[1]*$(syms_s[k]))
        for i = 2:k
            ex_line_s = :($ex_line_s + $i*a[$i]*$(syms_c[(k - i + 1)]))
            ex_line_c = :($ex_line_c - $i*a[$i]*$(syms_s[(k - i + 1)]))
        end
        ex_line_s = :($(syms_s[k + 1]) = ($ex_line_s)/$k)
        ex_line_c = :($(syms_c[k + 1]) = ($ex_line_c)/$k)
        ex_calc.args[2*k + 1] = ex_line_s
        ex_calc.args[2*k + 2] = ex_line_c
    end

    exout_s = :(($(syms_s[1]),))
    for i = 2:N
        push!(exout_s.args, syms_s[i])
    end

    exout_c = :(($(syms_c[1]),))
    for i = 2:N
        push!(exout_c.args, syms_c[i])
    end

    return quote
               Base.@_inline_meta
               $ex_calc
               return STaylor1{N,T}($exout_s), STaylor1{N,T}($exout_c)
            end
end

# Functions for STaylor1
@generated function tan(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:(N+1)])
    syms = Symbol[Symbol("c$i") for i in 1:N]
    syms2 = Symbol[Symbol("c2$i") for i in 1:N]

    ex_calc.args[1] = :($(syms[1]) = tan(a[0]))
    ex_calc.args[2] = :($(syms2[1]) = ($(syms[1]))^2)

    for k = 1:(N - 1)
        kT = convert(T, k)
        ex_line = :($(kT - 1)*a[$(k - 1)]*$(syms2[1]))
        @inbounds for i = 1:(k - 1)
            ex_line = :($ex_line + $(kT - i) * a[$(k - i)] * $(syms2[i + 1]))
        end
        ex_line = :(($ex_line)/$kT)
        ex_line = :($(syms[k + 1]) = $ex_line)
        ex_calc.args[k + 2] = ex_line
        #c2 = sqr(c)...
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

# Functions for STaylor1
@generated function asin(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(sym) = asin(a[0]))
    ex_calc.args[1] = ex_line

    #=
    for k in 1:(N-1)
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
    =#

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

# Functions for STaylor1
@generated function acos(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(sym) = acos(a[0]))
    ex_calc.args[1] = ex_line

    #=
    for k in 1:(N-1)
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
    =#

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

# Functions for STaylor1
@generated function atan(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(sym) = atan(a[0]))
    ex_calc.args[1] = ex_line

    #=
    for k in 1:(N-1)
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
    =#

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

sinh(a::STaylor1{N,T}) where {N, T <: Number} = sinhcosh(a)[1]
cosh(a::STaylor1{N,T}) where {N, T <: Number} = sinhcosh(a)[2]
@generated function sinhcosh(a::STaylor1{N,T}) where {N, T <: Number}

    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:(2*N)])

    syms_s = Symbol[Symbol("c$i") for i in 1:N]
    syms_c = Symbol[Symbol("c2$i") for i in 1:N]

    ex_line_s = :($(syms_s[1]) = sinh(a[0]))
    ex_line_c = :($(syms_c[1]) = cosh(a[0]))
    ex_calc.args[1] = ex_line_s
    ex_calc.args[2] = ex_line_c

    for k = 1:(N - 1)
        ex_line_s = :(a[1]*$(syms_c[k]))
        ex_line_c = :(a[1]*$(syms_s[k]))
        for i = 2:k
            ex_line_s = :($ex_line_s + $i*a[$i]*$(syms_c[(k - i + 1)]))
            ex_line_c = :($ex_line_c + $i*a[$i]*$(syms_s[(k - i + 1)]))
        end
        ex_line_s = :($(syms_s[k + 1]) = ($ex_line_s)/$k)
        ex_line_c = :($(syms_c[k + 1]) = ($ex_line_c)/$k)
        ex_calc.args[2*k + 1] = ex_line_s
        ex_calc.args[2*k + 2] = ex_line_c
    end

    exout_s = :(($(syms_s[1]),))
    for i = 2:N
        push!(exout_s.args, syms_s[i])
    end

    exout_c = :(($(syms_c[1]),))
    for i = 2:N
        push!(exout_c.args, syms_c[i])
    end

    return quote
               Base.@_inline_meta
               $ex_calc
               return STaylor1{N,T}($exout_s), STaylor1{N,T}($exout_c)
            end
end

# Functions for STaylor1
@generated function tanh(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(sym) = tanh(a[0]))
    ex_calc.args[1] = ex_line

    #=
    for k in 1:(N-1)
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
    =#

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
