function classify_sccs_by_cond(g,
    sccs=strongly_connected_components(g),
    cond=condensation(g, sccs)
)
    classes = []
    for (scc_i, scc) in enumerate(sccs)
        numstates = length(scc)
        ons = outneighbors(cond, scc_i)
        ins = inneighbors(cond, scc_i)

        if isempty(ons) # is an attracting components
            if isempty(ins) # is isolated
                push!(classes, numstates == 1 ? :tac : :iac)
            else
                push!(classes, numstates == 1 ? :sac : :gac)
            end
        elseif isempty(ins) # is a Garden of Eden
            push!(classes, numstates == 1 ? :goe1 : :goe)
        elseif length(ons) == length(ins) == 1
            push!(classes, numstates == 1 ? :spt1 : :spt)
        else
            push!(classes, numstates == 1 ? :gpt1 : :gpt)
        end
    end
    classes
end
classify_sccs_by_cond(gm::AbstractGraphModel, args...) = classify_sccs_by_cond(graph(gm), args...)
export classify_sccs_by_cond

function classify_acs_by_shape(g, acs=attracting_components(g))
    classes = []

    for ac in acs
        sg = induced_subgraph(g, ac)[1]
        iseq = true
        for e in edges(sg)
            if !has_edge(sg, e.dst, e.src)
                iseq = false
                break
            end
        end
        push!(classes, (iseq,))
    end

    classes
end
classify_acs_by_shape(gm::AbstractGraphModel, args...) = classify_acs_by_shape(graph(gm), args...)
export classify_acs_by_shape
