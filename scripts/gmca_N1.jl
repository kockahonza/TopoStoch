using DrWatson
@quickactivate "TopoStochSim"

using Revise

includet(scriptsdir("gmca_v2_5.jl"))

function sscur(; dofull=false)
    # global ca, sca, tm, rt, tm2, rt2, es, ssv, ssn, cs, sscv, sscn, sscv2, sscn2

    ca, sca = make_v2_5(1, 1; simplified=:both)
    tm, rt = etransmat_safe(ca; mat=true)
    tm2, rt2 = etransmat_safe(sca; mat=true)
    @assert isequal(tm, tm2)

    @variables a, b, c, d, e, f, g, h
    es = w_eigen(tm; safe=false)
    @assert es.evals[1] == 0
    ssv = es.evecs[1] * (a * d * e + a * d * f + b * c * f + b * d * f)
    ssn = w_simplify(sum(ssv); safe=false)

    cs = simplify.(calc_currents(tm, ssv); expand=true)
    sscv = w_simplify(ssubstitute(cs[3, 1], rt))
    if dofull
        sscn = w_simplify(ssubstitute(ssn, rt))
    else
        sscn = nothing
    end

    sscv2 = w_simplify(ssubstitute(cs[3, 1], rt2))
    sscn2 = w_simplify(ssubstitute(ssn, rt2))

    (; ca, sca, tm, rt, tm2, rt2, es, ssv, ssn, cs, sscv, sscn, sscv2, sscn2)
end

function eqss((; ssv, ssn, rt2))
    essv = w_simplify(ssubstitute(w_simplify(ssubstitute(ssv, rt2); full=true), get_eq_subs1()); full=true)
    essn = w_simplify(ssubstitute(w_simplify(ssubstitute(ssn, rt2); full=true), get_eq_subs1()); full=true)

    ess = w_simplify(essv / essn; full=true)

    (; essv, essn, ess)
end
eqss(_::Nothing) = eqss(sscur())

