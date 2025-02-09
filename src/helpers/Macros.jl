if debugging
    macro assert(args...)
        :($Base.@assert($(esc.(args)...)))
    end

    macro inbounds(arg)
        esc(arg)
    end
else
    macro assert(args...)
    end

    macro inbounds(arg)
        :($Base.@inbounds($(esc(arg))))
    end
end

macro verbose_info(str...)
    quote
        if $(esc(:verbose))
            println($(esc.(str)...))
            flush(stdout)
        end
    end
end

# Simpler version of FastClosure's @closure: We just require that the variable be marked at least once with an interpolation
# sign. In this way, there's precise control of what will be captured. The code here is almost identical to
# Base._lift_one_interp!, but we descend into macros (allowing for defining an escape list), and we also take care of duplicate
# interpolations. Once a variable was marked with an interpolation sign at any position, all occurrences are captured.
function _lift_interps!(e, escape_macros)
    letargs = Set{Any}()  # store the new gensymed arguments
    _lift_interps_helper(e, false, letargs, escape_macros) # Start out _not_ in a quote context (false)
    return letargs
end
_lift_interps_helper(v, _, _, escape_macros) = v
function _lift_interps_helper(expr::Expr, in_quote_context, letargs, escape_macros)
    if expr.head === :$
        if in_quote_context  # This $ is simply interpolating out of the quote
            # Now, we're out of the quote, so any _further_ $ is ours.
            in_quote_context = false
        else
            push!(letargs, :($(esc(expr.args[1])) = $(esc(expr.args[1]))))
            return expr.args[1] # Don't recurse into the lifted $() exprs
        end
    elseif expr.head === :quote
        in_quote_context = true   # Don't try to lift $ directly out of quotes
    elseif expr.head === :macrocall && expr.args[1] ∈ escape_macros
        return expr # Don't recur into escaped macro calls, since some other macros use $
    end
    for (i,e) in enumerate(expr.args)
        expr.args[i] = _lift_interps_helper(e, in_quote_context, letargs, escape_macros)
    end
    return expr
end

macro capture(firstarg, secondarg=nothing)
    if isnothing(secondarg)
        expr = firstarg
        escape_macros = Set{Symbol}()
    else
        (firstarg.head === :(=) && length(firstarg.args) === 2) || throw(ArgumentError("Invalid use of @capture"))
        firstarg.args[1] === :escape || throw(ArgumentError("Invalid keyword for @capture: $(firstarg.args[1])"))
        firstarg.args[2] isa Symbol || firstargs.args[2] ∈ (:vect, :tuple) ||
            throw(ArgumentError("Invalid value for keyword escape"))
        expr = secondarg
        escape_macros = Set{Symbol}(firstarg.args[2])
    end
    letargs = _lift_interps!(expr, escape_macros)
    quote
        let $(letargs...)
            $(esc(expr))
        end
    end
end

accumulate_vars!(v::Vector{Symbol}, e::Symbol) = push!(v, e)
function accumulate_vars!(v::Vector{Symbol}, e::Expr)
    e.head === :tuple || error("Variable assignments must be nested in a tuple-like syntax")
    for x in e.args
        accumulate_vars!(v, x)
    end
    return v
end

macro unroll(loop)
    loop.head === :for || error("@unroll requires a for loop")
    length(loop.args) == 2 || error("@unroll requires a simple for loop")
    vars, body = loop.args
    (vars.head === :(=) && length(vars.args) == 2) || error("Invalid for syntax")
    varnames, varvalues = vars.args
    (varvalues isa Expr && varvalues.head === :tuple) || error("@unroll requires tuple values")
    result = Expr(:block)
    sizehint!(result.args, length(varvalues.args))
    if varnames isa Symbol
        for varvalue in varvalues.args
            push!(result.args, :(let $varnames=$varvalue; $body end))
        end
    else
        allnames = Symbol[]
        accumulate_vars!(allnames, varnames)
        for varvalue in varvalues.args
            push!(result.args, :(let $(allnames...); $varnames = $varvalue; $body end))
        end
    end
    esc(result)
end