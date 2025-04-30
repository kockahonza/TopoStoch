using DrWatson
@quickactivate "TopoStochSim"

using GraphModels
using NonEqDigits

using DataFrames


"""
Classifies an attracting component as :single (1 vertex), :eq (all
edges bidirectional), or :neq (contains directional edges)
"""
function ac_eq_class(g, ac)
    if length(ac) == 1
        :single
    else
        ac_sg = induced_subgraph(g, ac)[1]
        for e in edges(ac_sg)
            if !has_edge(ac_sg, e.dst, e.src)
                return :neq
            end
        end
        :eq
    end
end
ac_eq_class(gm::AbstractGraphModel, ac) = ac_eq_class(graph(gm), ac)

"""
Checks if an attracting component forms a loop topology, ignoring edge directions.
"""
function ac_is_loop(g, ac)
    if length(ac) <= 1
        return false
    end

    ac_sg = induced_subgraph(g, ac)[1]
    for v in vertices(ac_sg)
        # Count total connections (in + out) for each vertex
        # In a loop, each vertex should have exactly two connections
        if length(all_neighbors(ac_sg, v)) != 2
            return false
        end
    end

    true
end
ac_is_loop(gm::AbstractGraphModel, ac) = ac_is_loop(graph(gm), ac)

"""
Categorizes attracting components into single vertices, equilibrium (eq) components,
non-equilibrium (neq) components, and their cyclic subcategories.
"""
function classify_acs_s1_stats(g, acs=attracting_components(g))
    acs_single = []
    acs_eq = []
    acs_neq = []
    acs_eqcyc = []
    acs_neqcyc = []

    for (ac_i, ac) in enumerate(acs)
        eq_class = ac_eq_class(g, ac)
        if eq_class == :single
            push!(acs_single, ac_i)
        else
            is_loop = ac_is_loop(g, ac)

            if eq_class == :eq
                push!(acs_eq, ac_i)
                if is_loop
                    push!(acs_eqcyc, ac_i)
                end
            elseif eq_class == :neq
                push!(acs_neq, ac_i)
                if is_loop
                    push!(acs_neqcyc, ac_i)
                end
            else
                throw(ErrorException("something went wrong in classify_acs_s1_stats"))
            end
        end
    end

    acs_single, acs_eq, acs_neq, acs_eqcyc, acs_neqcyc
end
classify_acs_s1_stats(gm::AbstractGraphModel, args...) = classify_acs_s1_stats(graph(gm), args...)

function make_s1_stats_df(Ns=3:6, codes=ca_ucodes_bydeg())
    df = DataFrame(;
        code=Int[],
        N=Int[],
        numstates=Int[],
        phi=Float64[],
        Nac=Int[],
        Nsingle=Int[],
        Neq=Int[],
        Nneq=Int[],
        Neqcyc=Int[],
        Nneqcyc=Int[],
        Mac=Int[],
        Msingle=Int[],
        Meq=Int[],
        Mneq=Int[],
        Meqcyc=Int[],
        Mneqcyc=Int[],
        Kac=Int[],
        Ksingle=Int[],
        Keq=Int[],
        Kneq=Int[],
        Keqcyc=Int[],
        Kneqcyc=Int[],
    )

    for code in codes
        for N in Ns
            ned = make_ca_ned(N, code; show=false)

            g = graph(ned)
            acs = attracting_components(g)

            ac_is_by_class = classify_acs_s1_stats(g, acs)
            ac_is_by_class = (1:length(acs), ac_is_by_class...)

            Ns_ = length.(ac_is_by_class)
            Ms = [maximum(ac_i -> length(acs[ac_i]), xx; init=0) for xx in ac_is_by_class]
            Ks = [sum(ac_i -> length(acs[ac_i]), xx; init=0) for xx in ac_is_by_class]

            numstates = 2^N
            phi = Ks[1] / numstates

            push!(df, (code, N, numstates, phi, Ns_..., Ms..., Ks...))
        end
    end

    df
end
