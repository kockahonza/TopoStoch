module NonEqDigits

using Reexport
@reexport using GraphModels

using DataFrames, CSV

# These are intended to be extended
import GraphModels: graph, numstates, allstates
import GraphModels: p_named_layouts, plotgm_layout_default, plotgm_kwarg_defaults
import Base: copy, broadcastable, show, display

################################################################################
# Base setup
################################################################################
abstract type Symmetry end
struct Chain <: Symmetry end
struct Loop <: Symmetry end
broadcastable(s::Symmetry) = Ref(s)
export Chain, Loop

"""
S is self explanatory, D is in essence the base aka the number of states per
digit, L is the length of the digit strings and F is the type of transition
weight values, mostly some float type. Safe is a Bool flag which is in the
type for performance.

This isn't captured in the type itself (for speed and sanity reasons), but
the states of this are along the lines of Vector{Int} with lenght L and each
Int being between >=0 and <D.
"""
struct NonEqDigitsGM{S<:Symmetry,D,L,F,Safe} <: AbstractGraphModel{F}
    graph::SimpleWeightedDiGraph{Int,F}

    kT::F # I keep this here so that it is consistent and usable throughout

    function NonEqDigitsGM(D, L, graph=nothing;
        symmetry=Loop,
        numtype=Float64,
        kT=1.0,
        safe=true
    )
        if isnothing(graph)
            graph = SimpleWeightedDiGraph{Int,numtype}(numstates(D, L))
        elseif nv(graph) != numstates(D, L)
            throw(ArgumentError(@sprintf "the number of vertices in the passed graph is not valid for D=%d and L=%d" D L))
        end
        if symmetry != Loop
            @warn "making a NonEqDigitsGM with symmetry not Loop, this is not fully supported yet!"
        end
        new{symmetry,D,L,numtype,safe}(graph, kT)
    end
end
graph(ned::NonEqDigitsGM) = ned.graph
export NonEqDigitsGM

################################################################################
# dealing with what the states are and finishing the interface implementation
################################################################################
numstates(D, L) = D^L
numstates(_::NonEqDigitsGM{S,D,L}) where {S,D,L} = numstates(D, L)

itostate_unsafe(i, D, L) = digits(i - 1; base=D, pad=L)
function itostate_safe(i, D, L)
    dstring = itostate_unsafe(i, D, L)
    if length(dstring) == L
        dstring
    else
        throw(ArgumentError(@sprintf "i=%d is not valid for D=%d and L=%d" i D L))
    end
end
function itostate(i, _::NonEqDigitsGM{S,D,L,F,Safe}) where {S,D,L,F,Safe}
    if Safe
        itostate_safe(i, D, L)
    else
        itostate_unsafe(i, D, L)
    end
end
export itostate

function statetoi_unsafe(dstring, D)
    i = 1
    di = 0
    for d in dstring
        i += d * D^di
        di += 1
    end
    i
end
function statetoi_safe(dstring, D, L)
    if (length(dstring) == L) && all(x -> 0 <= x < D, dstring)
        statetoi_unsafe(dstring, D)
    else
        throw(ArgumentError(@sprintf "dstring=%s is not valid for D=%d and L=%d" string(dstring) D L))
    end
end
function statetoi(st, _::NonEqDigitsGM{S,D,L,F,Safe}) where {S,D,L,F,Safe}
    if Safe
        statetoi_safe(st, D, L)
    else
        statetoi_unsafe(st, D)
    end
end
export statetoi

allstates_ned(D, L) = itostate_unsafe.(1:numstates(D, L), D, L)
allstates(_::NonEqDigitsGM{S,D,L}) where {S,D,L} = allstates_ned(D, L)

################################################################################
# Other utility operations
################################################################################
"""
given a digit index (assumed between 1 and L) returns a tuple of the left
and right neighboring indices, nothing if such a neighbor doesn't exist
"""
function neighbor_indices(i, _::NonEqDigitsGM{Loop,D,L}) where {D,L}
    if L == 1
        (nothing, nothing)
    elseif L == 2
        if i == 1
            (2, 2)
        else
            (1, 1)
        end
    else
        if i == 1
            (L, 2)
        elseif i == L
            (L - 1, 1)
        else
            (i - 1, i + 1)
        end
    end
end
export neighbor_indices

################################################################################
# Generic adding of transitions by mechanisms
################################################################################
"""
Add any general mechanism for the binary (D=2) class which only takes into
account up to NN interactions. Currently only supported for Loop.
"""
function add_mech_BNN!(ned::NonEqDigitsGM{S,2,L}, mu, K) where {S,L}
    for st_i in 1:numstates(ned)
        st = itostate(st_i, ned)
        for di in 1:L
            lni, rni = neighbor_indices(di, ned)
            lnd = st[lni] # left and right neighbor digits
            rnd = st[rni]
            if !isnothing(lnd) && !isnothing(rnd) # this might be static!
                kval = K[lnd+1, rnd+1]

                # this is done in this particular way so that it works with mu = +-Inf
                sigmoid_preval = 1.0 / (1.0 + exp(mu / ned.kT))
                if st[di] == 0
                    sigmoid_val = 1 - sigmoid_preval
                    new_d = 1
                else # implied that it is 1 from D==2
                    sigmoid_val = sigmoid_preval
                    new_d = 0
                end

                dest_st = copy(st)
                dest_st[di] = new_d

                inc_edge!(graph(ned), st_i, statetoi(dest_st, ned), kval * sigmoid_val)
            end
        end
    end
end
export add_mech_BNN!

neg_int(x::Int) = iszero(x) ? 1 : 0

# Should probably be a separate package but I can't be bothered
include("MolecularAutomata.jl")

################################################################################
# Adding other bits
################################################################################
include("plots_and_io.jl")
include("analysis.jl")
include("cellular_automata.jl")

end # module NonEqDigits
