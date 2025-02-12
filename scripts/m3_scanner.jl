using DrWatson
if isnothing(projectname())
    @quickactivate "TopoStochSim"
end

using JLD2
using StatsBase

using ComplexAllostery
using ComplexAllostery.Model3

using LinearAlgebra
using Base.Threads

################################################################################
# General and fast, threaded scan function for the simplified model 3
################################################################################
function _general_m3_scan_fast_bit(cca_fact, observables, data, do_ss, cis, r1s, r2s, r3s, alphas, ebs, cPs, cRs)
    # By default this is 11 but it seems much much better to parallelize here
    # instead of having a parallel eigensolver
    prev_blas_num_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    @threads for ci in cis
        r1_i, r2_i, r3_i, alpha_i, eb_i, cP_i, cR_i = ci.I
        cca = cca_fact(r1s[r1_i], r2s[r2_i], r3s[r3_i], alphas[alpha_i], ebs[eb_i], cPs[cP_i], cRs[cR_i])

        for (darray, observable) in zip(data, observables)
            if do_ss
                ss = supersteadystate(cca)
                darray[ci] = observable(cca, ss)
            else
                darray[ci] = observable(cca)
            end
        end
    end

    BLAS.set_num_threads(prev_blas_num_threads)
end
function general_m3_scan(N, observables, save_fname=nothing,
    obs_names=[f"obs_{i}" for i in 1:length(observables)];
    r1s=[1.0],
    r2s=[1.0],
    r3s=[1.0],
    alphas=[0.0],
    ebs=[0.0],
    cPs=[1.0],
    cRs=[1.0],
    do_ss=true,
    save_metadata=(;)
)
    ca = makeCAGM(N, vars_simplified())
    cca_fact::Function = make_factory(
        ca,
        [symvars.sim_rs; symvars.energies[3];
            symvars.concentrations[1]; symvars.redconcentration]
    )

    cis = CartesianIndices(length.((r1s, r2s, r3s, alphas, ebs, cPs, cRs)))

    # prep the data arrays
    fake_cca = cca_fact(1, 1, 1, 0, 0, 1, 1)
    if do_ss
        ss = supersteadystate(fake_cca)
    end
    data_types = []
    for observable in observables
        obs_return = if do_ss
            observable(fake_cca, ss)
        else
            observable(fake_cca)
        end
        push!(data_types, typeof(obs_return))
    end
    # kept as a tuple so that _general_m3_scan_fast_bit hopefully knows the types
    data = Tuple(Array{dt,length(cis.indices)}(undef, size(cis)) for dt in data_types)

    # run, separated to have a function barrier
    _general_m3_scan_fast_bit(cca_fact, observables, data, do_ss, cis, r1s, r2s, r3s, alphas, ebs, cPs, cRs)

    if !isnothing(save_fname)
        range_names = ["r1s", "r2s", "r3s", "alphas", "ebs", "cPs", "cRs"]
        obs_kwargs = Symbol.(obs_names) .=> data
        jldsave(save_fname; ca, r1s, r2s, r3s, alphas, ebs, cPs, cRs,
            range_names, observable_names=obs_names,
            obs_kwargs..., save_metadata...
        )
    end

    data
end

################################################################################
# Setting up different observables
################################################################################
"""
Constructs an observable for `general_m3_scan` that is just the ss probability
sum over sum nodes.
"""
struct ObsTotalAreaProbability
    area_indices::Vector{Int}
end
function (tap::ObsTotalAreaProbability)(cca, ss)
    sum(@view ss[tap.area_indices])
end

# Some really basic observables
function make_allligs_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x == B, st.occupations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end
function make_noligs_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x == 0, st.occupations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end
function make_alltense_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x == 1, st.conformations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end
function make_notense_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x != 1, st.conformations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end
function make_onboundary_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x == B, st.occupations)
            push!(ai, i)
        elseif all(x -> x == 0, st.occupations)
            push!(ai, i)
        elseif all(x -> x == 1, st.conformations)
            push!(ai, i)
        elseif all(x -> x != 1, st.conformations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end

# Going around the loop, very naive and simple measure of the loop
abstract type GRTLReturnType end
struct Matrix <: GRTLReturnType end
struct Mean <: GRTLReturnType end
struct MatrixAndMean <: GRTLReturnType end
struct List <: GRTLReturnType end

struct ObsGoingRoundTheLoop{RT<:GRTLReturnType,SS}
    start_nodes::Vector{Int}
    forwards::Vector{Vector{Int}}
    backwards::Vector{Vector{Int}}
    macro_positions::Vector{NTuple{2,Int}}
    function ObsGoingRoundTheLoop{RT,SS}(N, B=1) where {RT,SS}
        if !isa(SS, Bool)
            throw(ArgumentError("SS must be a Bool"))
        end
        start_nodes = []
        forwards = []
        backwards = []
        macro_positions = []
        for (st_i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, 2, B))
            numtense = calc_numofconf(st)
            numlig = calc_numligands(st)
            if (numlig == 0) || (numlig == N * B) || (numtense == 0) || (numtense == N)
                forward = []
                backward = []
                cst = copy(st)
                if (numlig == 0) && (numtense == 0)
                    for i in 1:N
                        cst.conformations[i] -= 1
                        push!(forward, statetoi(cst, N, 2, B))
                        cst.conformations[i] += 1

                        cst.occupations[i] += 1
                        push!(backward, statetoi(cst, N, 2, B))
                        cst.occupations[i] -= 1
                    end
                elseif (numlig == 0) && (numtense == N)
                    for i in 1:N
                        cst.occupations[i] += 1
                        push!(forward, statetoi(cst, N, 2, B))
                        cst.occupations[i] -= 1

                        cst.conformations[i] += 1
                        push!(backward, statetoi(cst, N, 2, B))
                        cst.conformations[i] -= 1
                    end
                elseif (numlig == N * B) && (numtense == N)
                    for i in 1:N
                        cst.conformations[i] += 1
                        push!(forward, statetoi(cst, N, 2, B))
                        cst.conformations[i] -= 1

                        cst.occupations[i] -= 1
                        push!(backward, statetoi(cst, N, 2, B))
                        cst.occupations[i] += 1
                    end
                elseif (numlig == N * B) && (numtense == 0)
                    for i in 1:N
                        cst.occupations[i] -= 1
                        push!(forward, statetoi(cst, N, 2, B))
                        cst.occupations[i] += 1

                        cst.conformations[i] -= 1
                        push!(backward, statetoi(cst, N, 2, B))
                        cst.conformations[i] += 1
                    end
                elseif numlig == 0
                    for i in 1:N
                        if st.conformations[i] != 1
                            cst.conformations[i] -= 1
                            push!(forward, statetoi(cst, N, 2, B))
                            cst.conformations[i] += 1
                        end

                        if st.conformations[i] == 1
                            cst.conformations[i] += 1
                            push!(backward, statetoi(cst, N, 2, B))
                            cst.conformations[i] -= 1
                        end
                    end
                elseif numlig == N * B
                    for i in 1:N
                        if st.conformations[i] == 1
                            cst.conformations[i] += 1
                            push!(forward, statetoi(cst, N, 2, B))
                            cst.conformations[i] -= 1
                        end

                        if st.conformations[i] != 1
                            cst.conformations[i] -= 1
                            push!(backward, statetoi(cst, N, 2, B))
                            cst.conformations[i] += 1
                        end
                    end
                elseif numtense == 0
                    for i in 1:N
                        if st.occupations[i] != 0
                            cst.occupations[i] -= 1
                            push!(forward, statetoi(cst, N, 2, B))
                            cst.occupations[i] += 1
                        end

                        if st.occupations[i] != B
                            cst.occupations[i] += 1
                            push!(backward, statetoi(cst, N, 2, B))
                            cst.occupations[i] -= 1
                        end
                    end
                elseif numtense == N
                    for i in 1:N
                        if st.occupations[i] != B
                            cst.occupations[i] += 1
                            push!(forward, statetoi(cst, N, 2, B))
                            cst.occupations[i] -= 1
                        end

                        if st.occupations[i] != 0
                            cst.occupations[i] -= 1
                            push!(backward, statetoi(cst, N, 2, B))
                            cst.occupations[i] += 1
                        end
                    end
                else
                    continue
                end
                push!(macro_positions, (numlig, numtense))
                push!(start_nodes, st_i)
                push!(forwards, forward)
                push!(backwards, backward)
            end
        end
        new{RT,SS}(start_nodes, forwards, backwards, macro_positions)
    end
end
function calc_rtl(rtl::ObsGoingRoundTheLoop{RT,false}, cca) where {RT}
    vals = fill(0.0, 1 + cca.N * cca.B, 1 + cca.N)
    for (st_i, fs, bs, pos) in zip(rtl.start_nodes, rtl.forwards, rtl.backwards, rtl.macro_positions)
        net_prop_rate = 0.0
        for f in fs
            net_prop_rate += get_weight(graph(cca), st_i, f)
        end
        for b in bs
            net_prop_rate -= get_weight(graph(cca), st_i, b)
        end
        net_prop_rate /= sum(on -> get_weight(graph(cca), st_i, on), outneighbors(graph(cca), st_i))
        vals[(pos .+ (1, 1))...] += net_prop_rate
    end
    vals
end
function calc_rtl(rtl::ObsGoingRoundTheLoop{RT,true}, cca, ss) where {RT}
    vals = fill(0.0, 1 + cca.N * cca.B, 1 + cca.N)
    for (st_i, fs, bs, pos) in zip(rtl.start_nodes, rtl.forwards, rtl.backwards, rtl.macro_positions)
        net_prop_rate = 0.0
        for f in fs
            net_prop_rate += get_weight(graph(cca), st_i, f)
        end
        for b in bs
            net_prop_rate -= get_weight(graph(cca), st_i, b)
        end
        net_prop_rate /= sum(on -> get_weight(graph(cca), st_i, on), outneighbors(graph(cca), st_i))
        vals[(pos .+ (1, 1))...] += ss[st_i] * net_prop_rate
    end
    vals
end
(rtl::ObsGoingRoundTheLoop{Matrix})(args...) = calc_rtl(rtl, args...)
(rtl::ObsGoingRoundTheLoop{Mean})(args...) = mean(calc_rtl(rtl, args...))
function (rtl::ObsGoingRoundTheLoop{MatrixAndMean})(args...)
    mat = calc_rtl(rtl, args...)
    mat, mean(mat)
end
function (rtl::ObsGoingRoundTheLoop{List})(cca)
    vals = Float64[]
    for (st_i, fs, bs) in zip(rtl.start_nodes, rtl.forwards, rtl.backwards)
        s = 0.0
        for f in fs
            s += get_weight(graph(cca), st_i, f)
        end
        for b in bs
            s -= get_weight(graph(cca), st_i, b)
        end
        s /= sum(get_weight.(Ref(graph(cca)), st_i, outneighbors(graph(cca), st_i)))
        push!(vals, s)
    end
    vals
end

function obs_ssavgnumlig(cca, ss)
    avgnumlig = 0.0
    for (st, p) in zip(allstates(cca), ss)
        avgnumlig += calc_numligands(st) * p
    end
    avgnumlig
end
function obs_ssavgnumofconf(cca, ss, args...)
    avgnumofconf = 0.0
    for (st, p) in zip(allstates(cca), ss)
        avgnumofconf += calc_numofconf(st, args...) * p
    end
    avgnumofconf
end
function obs_ssavgnumboundaries(cca, ss)
    avgnumboundaries = 0.0
    for (st, p) in zip(allstates(cca), ss)
        avgnumboundaries += calc_numboundaries(st, cca) * p
    end
    avgnumboundaries
end

################################################################################
# Local testing
################################################################################
function scantest1()
    observables = [make_onboundary_obs(4)]

    general_m3_scan(4, observables;
        r1s=10 .^ LinRange(-5, 3, 15),
        # r2s=10 .^ LinRange(-5, 3, 15),
        alphas=[0, 0.5],
        ebs=LinRange(0, 6, 5),
    )
end

function rtltest(N=4)
    observables = [ObsGoingRoundTheLoop{Matrix,false}(N), ObsGoingRoundTheLoop{Mean,false}(N)]

    general_m3_scan(N, observables;
        r1s=10 .^ LinRange(-5, 3, 15),
        # r2s=10 .^ LinRange(-5, 3, 15),
        alphas=[0, 0.5],
        ebs=LinRange(0, 6, 5),
        do_ss=false
    )
end
