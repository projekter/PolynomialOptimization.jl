function PolynomialOptimization.tighten_minimize_l1(::Val{:Mosek}, spmat, rhs)
    maketask() do task
        len = size(spmat, 2)
        # For the ℓ₁ norm, we need the absolute values, so another len variables.
        appendvars(task, 2len)
        putvarboundsliceconst(task, 1, 2len +1, MSK_BK_FR, -Inf, Inf)
        # Minimize the absolute values
        putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE)
        putclist(task, collect(len+1:2len), fill(1., len))
        # Make sure the original variables satisfy spmat*x = rhs
        lensp = size(spmat, 1)
        appendcons(task, lensp + 2len)
        putacolslice(task, 1, len +1, spmat)
        putconboundslice(task, 1, lensp +1, fill(MSK_BK_FX, lensp), rhs, rhs)
        # And that the others are larger or equal to the absolute values.
        for (j, i) in enumerate(lensp+1:2:lensp+2len)
            putarow(task, i, [j, j + len], [1., -1.])
            putarow(task, i +1, [j, j + len], [-1., -1.])
        end
        putconboundsliceconst(task, lensp +1, lensp + 2len +1, MSK_BK_UP, -Inf, 0.)
        optimize(task)
        status = getsolsta(task, MSK_SOL_BAS)
        if status === MSK_SOL_STA_OPTIMAL
            return getxxslice(task, MSK_SOL_BAS, 1, len +1)
        else
            throw(SingularException(0))
        end
    end
end