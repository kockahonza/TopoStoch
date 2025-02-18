module NonEqDigits

using Reexport
@reexport using GraphModels

using Printf

# These are intended to be extended
import GraphModels: graph, numstates, allstates
import GraphModels: p_named_layouts, p_do_layout, plotgm_kwarg_defaults
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
weight values, mostly some float type.
"""
struct NonEqDigitsGM{S<:Symmetry,D,L,Safe,F} <: AbstractGraphModel{F}

end

function dstringtoi(dstring, D, L)
    i = 1
    di = 0
    for d in dstring
        i += d * D^di
        di += 1
    end
    i
end
export dstringtoi
struct NEDState{D,L}
    i::Int # the number this state corresponds to
    dstring::SVector{L,Int} # of lenght L and guaranteed to only contain digits 1 to D

    # constructing from a dstring
    function NEDState{D,L}(dstring) where {D,L}
        if all(x -> 0 <= x < D, dstring)
            new{D,L}(1, dstring)
        else
            throw(ArgumentError(@sprintf "cannot construct NEDState{%d,%d} with given digit dstring %s" D L string(dstring)))
        end
    end
    NEDState{D}(dstring) where {D} = NEDState{D,length(dstring)}(dstring)

    # constructing from i
end
export NEDState

numstates(D, L) = D^L
numstates(_::NonEqDigitsGM{S,D,L}) where {S,D,L} = numstates(D, L)

itostate_unsafe(i, D, L) = digits(i-1; base=D, pad=L)
function itostate(i, neq::NonEqDigitsGM{S,D,L,Safe}) where {S,D,L,Safe}
    if Safe
        if (i < 1) || (i > D^L - 1)
            throw(ArgumentError(@sprintf "i=%d is out of range for D=%d and L=%d" i D L))
        end
    end
    itostate_unsafe(i, D, L)
end

end # module NonEqDigits
