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

# Functions for STaylor1
@generated function sin(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(sym) = sin(a[0]))
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
@generated function cos(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(sym) = cos(a[0]))
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
@generated function tan(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(sym) = tan(a[0]))
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

# Functions for STaylor1
@generated function sinh(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(sym) = sinh(a[0]))
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
@generated function cosh(a::STaylor1{N,T}) where {N, T <: Number}
    ex_calc = quote end
    append!(ex_calc.args, Any[nothing for i in 1:N])
    syms = Symbol[Symbol("c$i") for i in 1:N]

    sym = syms[1]
    ex_line = :($(sym) = cosh(a[0]))
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
