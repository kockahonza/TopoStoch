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

"""
Returns the size/length of the only non-trivial (more than 1 state) attracting
components if possible. Returns 0 otherwise.
"""
function calc_onlynontrivaclen(g)
    onontrivac = get_onlynontrivac(g)
    if !isnothing(onontrivac)
        length(onontrivac)
    else
        0
    end
end
export get_onlynontrivac, calc_onlynontrivaclen
