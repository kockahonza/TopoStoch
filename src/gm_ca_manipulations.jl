################################################################################
# ca/graph/matrix symbolics manipulation using Wolfram
################################################################################
includet(srcdir("letter_codes.jl"))

function substitute_safenames(obj)
    terms = Dict{Num,Num}()
    code = Code()
    for var in get_variables(obj)
        terms[var] = Symbolics.variable(repr(code))
        increment!(code)
    end
    rev_terms = map(reverse, collect(terms))
    ssubstitute(obj, terms), rev_terms
end

function safemat_general(matrix)
    simple_matrix = similar(matrix)
    rev_terms = Dict{Num,Num}()
    code = Code()
    for i in eachindex(matrix)
        v = matrix[i]
        if iszero(v)
            simple_matrix[i] = v
        else
            var = Symbolics.variable(repr(code))
            simple_matrix[i] = var
            increment!(code)
            rev_terms[var] = v
        end
    end
    simple_matrix, rev_terms
end

function etransmat_safe(gm::AbstractGraphModel, args...; kwargs...)
    tm = transmat(gm, args...; kwargs...)
    setm, rt = safemat_general(tm)
    etransmat!(setm)
    setm, rt
end

function prep_for_wolfram(obj)
    sobj, rt = substitute_safenames(obj)
    SymbolicsMathLink.expr_to_mathematica(sobj), rt
end

function make_safe(safe, obj)
    if safe
        sobj, rt = substitute_safenames(obj)
        sobj, (srlst -> ssubstitute(srlst, rt))
    else
        obj, identity
    end
end

function w_simplify(obj; full=true, safe=true)
    obj, desafe = make_safe(safe, obj)
    cmd = full ? "FullSimplify" : "Simplify"
    rslt = wcall.(cmd, obj)
    desafe(rslt)
end

function w_eigen(matrix; safe=true)
    matrix, desafe = make_safe(safe, matrix)
    rslt = wcall("Eigensystem", collect(transpose(matrix)))
    (evals=desafe(rslt[1]), evecs=desafe.(rslt[2]))
end

function w_test_eigen(matrix; eigen_func=w_eigen, safe=false)
    (evals, evecs) = eigen_func(matrix; safe)
    residues = []
    for (eval, evec) in zip(evals, evecs)
        coeffs = (matrix * evec) ./ evec
        push!(residues, w_simplify(coeffs .- eval; safe))
    end
    residues
end

# Particular, less general functions
function ca_eigen1(ca::ComplexAllosteryGM; simplify=false)
    setm, rt = etransmat_safe(ca)
    esys = w_eigen(setm; safe=false)
    if simplify
        esys = map(x -> w_simplify(x; safe=false), esys)
    end
    esys, rt
end

function ca_eigen2(ca::ComplexAllosteryGM; simplify=false)
    etm = etransmat(ca)
    setm, rt = substitute_safenames(etm)
    esys = w_eigen(setm; safe=false)
    if simplify
        esys = map(x -> w_simplify(x; safe=false), esys)
    end
    esys, rt
end

# Shorthands
ssn = substitute_safenames
function towolfram(term::Num)
    SymbolicsMathLink.expr_to_mathematica(term)
end
function towolfram(term::AbstractMatrix)
    SymbolicsMathLink.expr_to_mathematica(collect(transpose(term)))
end
