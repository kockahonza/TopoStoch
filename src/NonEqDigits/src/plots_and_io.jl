function plotgm_kwarg_defaults(ned::NonEqDigitsGM{S,D,L}) where {S,D,L}
    rslt = Dict{Symbol,Any}()
    rslt[:fnlabels] = join
    rslt[:felabels] = false
    rslt[:n_ss_size] = false
    rslt[:ss] = false # until I fix supersteadystate for many sss

    if (D == 2) && (L == 3)
        rslt[:layout] = collect(allstates(ned))
    end

    rslt
end

function plotgm_kwarg_defaults(ma::MolAut)
    rslt = Dict{Symbol,Any}()
    rslt[:fnlabels] = join
    rslt[:felabels] = false
    rslt[:n_ss_size] = false
    rslt[:ss] = false # until I fix supersteadystate for many sss

    # if (ma.L == 3)
    #     rslt[:layout] = collect(labels(ma.mg))
    # end

    rslt
end
