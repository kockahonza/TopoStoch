module WSymbolics

using GraphModels

using Reexport
@reexport using SymbolicsMathLink

################################################################################
# Base
################################################################################
towolfram = SymbolicsMathLink.expr_to_mathematica

function prep_for_wolfram(obj)
    sobj, rt = substitute_safenames(obj)
    towolfram(sobj), rt
end
export prep_for_wolfram

function make_safe(safe, obj)
    if safe
        sobj, rt = substitute_safenames(obj)
        sobj, (srlst -> ssubstitute(srlst, rt))
    else
        obj, identity
    end
end

function w_simplify(obj; full=false, safe=true)
    obj, desafe = make_safe(safe, obj)
    cmd = full ? "FullSimplify" : "Simplify"
    rslt = wcall.(cmd, obj)
    desafe(rslt)
end
export w_simplify

function swcall(cmd, obj; safe=true)
    obj, desafe = make_safe(safe, obj)
    rslt = wcall.(cmd, obj)
    desafe(rslt)
end
export swcall

# TODO: This should be updated to use Eigensystem
function w_eigen(matrix; safe=true)
    matrix, desafe = make_safe(safe, matrix)
    rslt = wcall("Eigensystem", collect(matrix))
    (evals=desafe(rslt[1]), evecs=desafe.(rslt[2]))
end
export w_eigen

function w_test_eigen(matrix; eigen_func=w_eigen, safe=false)
    (evals, evecs) = eigen_func(matrix; safe)
    residues = []
    for (eval, evec) in zip(evals, evecs)
        coeffs = (matrix * evec) ./ evec
        push!(residues, w_simplify(coeffs .- eval; safe))
    end
    residues
end

################################################################################
# Potentially load other things using Requires
################################################################################
using Requires

function __init__()
    @require ComplexAllostery = "cb963ee1-309c-4455-b8c2-80df2e9615da" include("ext_ComplexAllostery.jl")
end

end # module WSymbolics
