"""
Gillespie algorithm that should actually be fast and parallelizeable, save
output in jld2 files which includes the whole trajectory or observables,
optionally does some plotting
"""
module GillespieSim

using GraphModels

using Base.Threads
using Dates
using Random
using StatsBase
using HDF5
using JLD2

import Base: close

################################################################################
# Main, performance-critical Gillespie struct AbstractGillespieSaver API
################################################################################
"""
T and F correspond to some integer and floating point number types.
Realistically I might have as well made them concrete Int and FLoat64 but
I guess this is a little more general.
"""
mutable struct Gillespie{T,F,R<:AbstractRNG}
    # General
    graph::SimpleWeightedDiGraph{T,F}
    # Current state
    cur_time::F
    cur_state::T
    # Internals to make things stable and fast
    rng::R
    r0s::Vector{F}
    cached_neighbors::Vector{Vector{T}}
    cached_neighbor_rates::Vector{Weights{F,F,Vector{F}}}
    function Gillespie(g::SimpleWeightedDiGraph{T,F}, start_time=0.0, start_state=1;
        rng=nothing
    ) where {T,F}
        if isnothing(rng)
            rng = Xoshiro()
        elseif isa(rng, Number)
            rng = Xoshiro(rng)
        elseif !isa(rng, AbstractRNG)
            throw(ArgumentError(f"kwarg rng must be nothing a number of an AbstractRNG but is a {typeof(rng)}"))
        end

        r0s = Vector{F}(undef, nv(g))
        cns = Vector{Vector{T}}(undef, nv(g))
        cws = Vector{Weights{F,F,Vector{F}}}(undef, nv(g))
        for i in 1:nv(g)
            r0s[i] = sum(get_weight.(Ref(g), i, neighbors(g, i)))
            cns[i] = neighbors(g, i)
            cws[i] = Weights(get_weight.(Ref(g), i, cns[i]))
        end

        new{T,F,typeof(rng)}(g, start_time, start_state, rng, r0s, cns, cws)
    end
end
Gillespie(gm::AbstractGraphModel, args...; kwargs...) = Gillespie(graph(gm), args...; kwargs...)
export Gillespie

"""
Some saver that will work with `Gillespie{T,F}` types. It has to at least save
a time series of type `F` and which can be accessed through `times`.
"""
abstract type AbstractGillespieSaver{T,F} end
function save!(saver::AbstractGillespieSaver{T,F}, _::Gillespie{T,F}) where {T,F}
    throw(ErrorException(f"No method of \"save!\" was provided for type \"{typeof(saver)}\""))
end
function close(saver::AbstractGillespieSaver)
    throw(ErrorException(f"No method of \"close\" was provided for type \"{typeof(saver)}\""))
end
function times(saver::AbstractGillespieSaver)
    throw(ErrorException(f"No method of \"times\" was provided for type \"{typeof(saver)}\""))
end
function last_saved_time(saver::AbstractGillespieSaver{T,F})::Union{Nothing,F} where {T,F}
    times_ = times(saver)
    if length(times_) > 1
        times_[end]
    else
        nothing
    end
end
function valid_saver(saver::AbstractGillespieSaver, gs::Gillespie)
    last_time = last_saved_time(saver)
    if isnothing(last_time) || isapprox(gs.cur_time, last_time)
        true
    else
        false
    end
end

################################################################################
# Concrete savers
################################################################################
struct SimpleH5Saver{T,F} <: AbstractGillespieSaver{T,F}
    destination::HDF5.H5DataStore
    time_dataset::HDF5.Dataset
    path_dataset::HDF5.Dataset

    function SimpleH5Saver{T,F}(destination) where {T,F}
        time_dataset = create_dataset(destination, "time", datatype(F), ((0,), (-1,)), chunk=(1,))
        path_dataset = create_dataset(destination, "path", datatype(T), ((0,), (-1,)), chunk=(1,))

        new{T,F}(destination, time_dataset, path_dataset)
    end
end
function SimpleH5Saver(_::Gillespie{T,F}, destination) where {T,F}
    SimpleH5Saver{T,F}(destination)
end
function SimpleH5Saver(_::SimpleWeightedDiGraph{T,F}, destination) where {T,F}
    SimpleH5Saver{T,F}(destination)
end
function SimpleH5Saver(gm::AbstractGraphModel, destination)
    SimpleH5Saver(graph(gm), destination)
end
export SimpleH5Saver
function save!(sgs::SimpleH5Saver{T,F}, gs::Gillespie{T,F}) where {T,F}
    new_length = length(sgs.time_dataset) + 1

    HDF5.set_extent_dims(sgs.time_dataset, (new_length,))
    sgs.time_dataset[new_length] = gs.cur_time

    HDF5.set_extent_dims(sgs.path_dataset, (new_length,))
    sgs.path_dataset[new_length] = gs.cur_state
end
function close(sgs::SimpleH5Saver)
    if isa(sgs.destination, HDF5.File)
        close(sgs.destination)
        true
    else
        false
    end
end
times(sgs::SimpleH5Saver) = sgs.time_dataset

################################################################################
# Running, contains the actual algorithm
################################################################################
"""
Returns the time until next transition and the destination
"""
function next_step(gs::Gillespie{T,F}) where {T,F}
    delta_time = randexp(gs.rng, F) / gs.r0s[gs.cur_state]
    ns = gs.cached_neighbors[gs.cur_state]
    if isempty(ns)
        nothing
    else
        next_state = sample(
            gs.rng,
            gs.cached_neighbors[gs.cur_state],
            gs.cached_neighbor_rates[gs.cur_state]
        )
        (delta_time, next_state)
    end
end

"""
Only used for verification of the above. Numerically seems to agree, so good!
"""
function next_step_direct(gs::Gillespie{T,F}) where {T,F}
    r1 = rand(gs.rng, F)
    r2 = rand(gs.rng, F)

    delta_time = -log(r1) / gs.r0s[gs.cur_state]
    next_state_i = 1
    while sum(
        i -> get_weight(gs.graph, gs.cur_state, i),
        gs.cached_neighbors[gs.cur_state][begin:next_state_i]
    ) <= r2 * gs.r0s[gs.cur_state]
        next_state_i += 1
    end

    (delta_time, gs.cached_neighbors[gs.cur_state][next_state_i])
end

function run_gillespie!(gs::Gillespie;
    saver=nothing,
    callback=nothing,
    logging=true,
    # Exit conditions
    exit_steps=nothing,
    exit_sim_deltatime=nothing,
    exit_sim_time=nothing,
    exit_real_time=nothing
)
    # Do some checks
    if !isnothing(saver)
        if !valid_saver(saver, gs)
            throw(ErrorException("The passed saver and Gillespie are not compatible, cannot start run"))
        end
    end
    if isnothing(exit_steps) && isnothing(exit_sim_deltatime) && isnothing(exit_sim_time) && isnothing(exit_real_time)
        throw(ArgumentError("cannot run `run_gillespie!` without at least one exit condition"))
    end

    # Prepare exit conditions
    sim_steps = 0
    sim_time = 0.0
    if !isnothing(exit_real_time)
        exit_real_time_ = if isa(exit_real_time, Number)
            Second(exit_real_time)
        else
            exit_real_time
        end
        start_real_time = now()
    end

    exit = false
    while !exit
        # Make a step
        ns = next_step(gs)
        if isnothing(ns)
            return
        end

        (delta_time, next_state) = ns
        gs.cur_time += delta_time
        gs.cur_state = next_state
        sim_steps += 1
        sim_time += delta_time

        # Do saving/callbacks
        if !isnothing(saver)
            save!(saver, gs)
        end

        if !isnothing(callback)
            callback(gs)
        end

        # Check exit conditions
        if !isnothing(exit_steps)
            if sim_steps >= exit_steps
                if logging
                    @info "Exited due to max simulation number of steps met"
                end
                exit = true
            end
        end
        if !isnothing(exit_sim_deltatime)
            if sim_time >= exit_sim_deltatime
                if logging
                    @info "Exited due to simulation timepoint reached time met"
                end
                exit = true
            end
        end
        if !isnothing(exit_sim_time)
            if gs.cur_time >= exit_sim_time
                if logging
                    @info "Exited due to max simulation time met"
                end
                exit = true
            end
        end
        if !isnothing(exit_real_time)
            if (now() - start_real_time) >= exit_real_time_
                if logging
                    @info "Exited due to max real world time met"
                end
                exit = true
            end
        end
    end
end
export run_gillespie!

################################################################################
# Particular usage cases which run gillespie sims and also do other things like
# save the GraphModel etc. along with loader functions. These should be used
# primarily.
################################################################################
function _default_gm_datadir(gm::AbstractGraphModel)
    datadir("scratch", f"{string(typeof(gm))}_{savename(gm)}")
end

# Running/loading a single sim
function run_single_gillespie(
    gm::AbstractGraphModel{<:AbstractFloat},
    start_node,
    save_path::String=_default_gm_datadir(gm);
    override=false,
    rng=nothing,
    kwargs...
)
    if ispath(save_path)
        if !override
            throw(ArgumentError(f"save_path {save_path} already exists"))
        else
            rm(save_path; recursive=true)
        end
    end
    mkpath(save_path)
    @info f"Starting gillespie single run at {save_path} starting from {start_node}"

    # Prep the gillespie run
    gs = Gillespie(gm, 0.0, start_node; rng)
    # Prep the saver
    filepath = joinpath(save_path, "gillespie_data.h5")
    saver = SimpleH5Saver(gm, h5open(filepath, "w"))
    # Save the gm
    save_object(joinpath(save_path, "gm.jld2"), gm)

    run_gillespie!(gs; saver, kwargs...)

    if !close(saver)
        @warn f"Could not close saver"
    end
end
export run_single_gillespie
function load_single_gillespie(save_path;
    gm_filename="gm.jld2",
    gillespie_filename="gillespie_data.h5"
)
    load_object(joinpath(save_path, gm_filename)), h5open(joinpath(save_path, gillespie_filename))
end
export load_single_gillespie

# Running/loading one sim starting from each node in gm
"""
Will run multiple Gillespie simulations, one starting at each possible state
all the runs are saved independently a newly created directory `save_path`
errors if path already exists.
"""
function run_full_gillespie_ensemblesim(
    gm::AbstractGraphModel{<:AbstractFloat},
    save_path::String=_default_gm_datadir(gm);
    override=false,
    rng=nothing,
    kwargs...
)
    if ispath(save_path)
        if !override
            throw(ArgumentError(f"save_path {save_path} already exists"))
        else
            rm(save_path; recursive=true)
        end
    end
    mkpath(save_path)
    @info f"Starting gillespie full ensemble run at {save_path}"

    # Save the gm
    save_object(joinpath(save_path, "gm.jld2"), gm)

    @threads for start_node in 1:numstates(gm)
        gs = Gillespie(gm, 0.0, start_node; rng)
        filepath = joinpath(save_path, f"start_{start_node}.h5")
        saver = SimpleH5Saver(gm, h5open(filepath, "w"))

        run_gillespie!(gs; saver, kwargs...)

        if !close(saver)
            @warn f"Could not close saver from start node {start_node}"
        end
        @info f"Finished running from start node {start_node}"
        flush(stdout)
        GC.gc()
    end
end
export run_full_gillespie_ensemblesim
function load_full_gillespie_ensemblesim(save_path, i=nothing;
    gm_filename="gm.jld2",
    save_prefix="start_"
)
    gm = try
        load_object(joinpath(save_path, gm_filename))
    catch
        nothing
    end
    loader_or_file = if isnothing(i)
        function (i::Integer)
            f = h5open(joinpath(save_path, save_prefix * string(i) * ".h5"))
            f["time"], f["path"]
        end
    else
        h5open(joinpath(save_path, save_prefix * string(i) * ".h5"))
    end
    gm, loader_or_file
end
export load_full_gillespie_ensemblesim

end
