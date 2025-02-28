"""
Returns the the only non-trivial (more than 1 state) strongly connected
component if possible and nothing otherwise.
"""
function get_onlynontrivscc(g)
    nontrivacs = filter(x -> length(x) != 1, strongly_connected_components(g))
    if length(nontrivacs) == 1
        nontrivacs[1]
    else
        nothing
    end
end
get_onlynontrivscc(gm::AbstractGraphModel{<:AbstractFloat}) = get_onlynontrivscc(graph(gm))
get_onlynontrivscclen(g) = length(something(get_onlynontrivscc(g), ()))
export get_onlynontrivscc, get_onlynontrivscclen

"""
Returns the the only non-trivial (more than 1 state) attracting component
if possible and nothing otherwise.
"""
function get_onlynontrivac(g)
    nontrivacs = filter(x -> length(x) != 1, attracting_components(g))
    if length(nontrivacs) == 1
        nontrivacs[1]
    else
        nothing
    end
end
get_onlynontrivac(gm::AbstractGraphModel{<:AbstractFloat}) = get_onlynontrivac(graph(gm))
get_onlynontrivaclen(g) = length(something(get_onlynontrivac(g), ()))
export get_onlynontrivac, get_onlynontrivaclen

"""
Finds all strongly connected components of gm that have only one element/state
and that are not connected to any other scc.
"""
function findsingularsccs(
    gm::AbstractGraphModel{<:AbstractFloat}, sccs::Vector,
    cond=condensation(graph(gm), sccs) # can be passed for speed
)
    singular_sccs = []
    nonsingular_sccs = []
    for scci in 1:length(sccs)
        if (length(sccs[scci]) == 1) && isempty(all_neighbors(cond, scci))
            push!(singular_sccs, scci)
        else
            push!(nonsingular_sccs, scci)
        end
    end
    singular_sccs, nonsingular_sccs
end
export findsingularsccs
