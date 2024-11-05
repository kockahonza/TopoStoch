################################################################################
# ca/graph manipulations
################################################################################
function filter_edges!(graph, threshold)
    for e in edges(graph)
        if weight(e) < threshold
            rem_edge!(graph, e)
        end
    end
end
function filter_edges!(ca::ComplexAllosteryGM{S,F}, args...) where {S,F}
    if F == Num
        throw(ArgumentError(f"function filter_edges! is only valid for ComplexAllosteryGM types with a concrete weight type - aka F != Num"))
    end
    filter_edges!(ca.graph, args...)
end

function keep_best_only!(graph, threshold=1e-10)
    for vertex in 1:nv(graph)
        neightbors_ = copy(neighbors(graph, vertex))

        max_rate = 0.0
        for n in neightbors_
            nrate = get_weight(graph, vertex, n)
            if nrate > max_rate
                max_rate = nrate
            end
        end

        to_remove = []
        for n in neightbors_
            nrate = get_weight(graph, vertex, n)
            if (max_rate - nrate) > threshold
                push!(to_remove, n)
            end
        end
        for n in to_remove
            rem_edge!(graph, vertex, n)
        end
    end
end
function keep_best_only!(ca::ComplexAllosteryGM{S,F}, args...) where {S,F}
    if F == Num
        throw(ArgumentError(f"function keep_best_only! is only valid for ComplexAllosteryGM types with a concrete weight type - aka F != Num"))
    end
    keep_best_only!(ca.graph, args...)
end

function copyand(f!, args...; kwargs...)
    function (obj)
        cobj = copy(obj)
        f!(cobj, args...; kwargs...)
        cobj
    end
end

################################################################################
# ca/graph/matrix symbolics manipulation using Wolfram
################################################################################
function substitute_safenames(obj)
    terms = Dict{Num,Num}()
    letter = 'a'
    for var in get_variables(obj)
        terms[var] = Symbolics.variable(letter)
        letter += 1
        if letter > 'z'
            throw(ErrorException("object has too many variables, ran out of alphabet"))
        end
    end
    rev_terms = map(reverse, collect(terms))
    ssubstitute(obj, terms), rev_terms
end

function substitute_safenames2(matrix)
    simple_matrix = similar(matrix)
    rev_terms = Dict{Num,Num}()
    letter = 'a'
    for i in eachindex(matrix)
        v = matrix[i]
        if iszero(v)
            simple_matrix[i] = v
        else
            var = Symbolics.variable(letter)
            simple_matrix[i] = var
            letter += 1
            rev_terms[var] = v
        end
    end
    simple_matrix, rev_terms
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

function w_eigen(ca::ComplexAllosteryGM; simplify=false)
    tm = transmat(ca)
    stm, rt = substitute_safenames2(tm)
    esys = w_eigen(stm; safe=false)
    if simplify
        esys = map(x -> w_simplify(x; safe=false), esys)
    end
    map(x -> ssubstitute(x, rt), esys)
end
