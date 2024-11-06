################################################################################
# Getting the transition matrices in a general way
################################################################################
"""
The standard Markov chain/Master equation transition matrix where the
ij th element corresponds to a transition/rate from j to i.
"""
function transmat(gm::AbstractGraphModel; mat=false)
    tm = transpose(adjacency_matrix(graph(gm)))
    mat ? Matrix(tm) : sparse(tm)
end

"""
Returns the modified transition matrix that can be used to get the time
derivative as matrix multiplication only. See latex_notes/MasterEq this
is W. Follows the same convention as `transmat`.
"""
function etransmat!(tm::AbstractMatrix)
    n = size(tm)[1]
    for i in 1:n
        tm[i, i] -= sum(tm[k, i] for k in 1:n)
    end
end
function etransmat(tm::AbstractMatrix)
    ctm = copy(tm)
    etransmat!(ctm)
    ctm
end
function etransmat(gm::AbstractGraphModel, args...; kwargs...)
    tm = transmat(gm, args...; kwargs...)
    etransmat!(tm)
    tm
end

################################################################################
# Doing stuff with the etransmat
################################################################################
function fixevec(evec::AbstractVector)
    if all(evec .<= 0)
        true, -evec
    elseif all(evec .>= 0)
        true, evec
    else
        false, evec
    end
end
function smarteigen(gm::AbstractGraphModel; norm=true, normthreshold=1e-8)
    esys = eigen(etransmat(gm; mat=true); sortby=abs)
    n = length(esys.values)

    evecs = Vector{Vector{Float64}}(undef, n)
    valid = BitVector(undef, n)
    for i in 1:n
        valid[i], evecs[i] = fixevec(esys.vectors[:, i])
        evecsum = sum(evecs[i])
        if evecsum < 0.0
            evecs[i] .*= -1.0
            evecsum *= -1.0
        end
        if norm && (evecsum > normthreshold)
            evecs[i] ./= evecsum
        end
    end

    (; evals=esys.values, evecs, valid)
end

################################################################################
# Editing/filtering the graph
################################################################################
function filter_edges!(graph::AbstractFloatGraph, threshold)
    for e in edges(graph)
        if weight(e) < threshold
            rem_edge!(graph, e)
        end
    end
end
function filter_edges!(ca::AbstractGraphModel, args...)
    filter_edges!(graph(ca), args...)
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
function keep_best_only!(ca::AbstractGraphModel, args...)
    keep_best_only!(graph(ca), args...)
end

# util func for the modifying functions above
function copyand(f!, args...; kwargs...)
    function (obj)
        cobj = copy(obj)
        f!(cobj, args...; kwargs...)
        cobj
    end
end


# TODO: FIX: This is not finished!!
function find_cycles(graph, next_vertex_func=next_vertex_choose_max)
    cycles = []
    cycle_map = Vector{Union{Nothing,Int}}(nothing, nv(graph))
    for vertex in 1:nv(graph)
        if isnothing(cycle_map[vertex])
            cycle = []
            true_cycle_start_i = findfirst(x -> x == vertex, cycle)
            while isnothing(true_cycle_start_i)
                push!(cycle, vertex)
                vertex = next_vertex_func(graph, vertex)
                true_cycle_start_i = findfirst(x -> x == vertex, cycle)
            end
            for v in cycle
                cycle_map[v] = length(cycles) + 1
            end
            cycle = cycle[true_cycle_start_i:end]
            push!(cycles, cycle)
        end
    end
    (cycles, cycle_map)
end
