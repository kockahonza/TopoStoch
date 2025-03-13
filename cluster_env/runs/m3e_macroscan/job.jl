using DrWatson
@quickactivate "cluster_env"

include("../../../scripts/m3_scanner.jl")

function macroscan(N=4; testrun=false)
    observables = [obs_ssavgnumlig, obs_ssavgenergy, obs_ssavgnumofconf, obs_ssavgnumboundaries]
    obs_names = ["numlig", "energy", "numtense", "numboundaries"]

    fname = f"./macroscan_N{N}.jld2"
    if ispath(fname)
        @error f"Path {fname} already exists, exiting so that it is not overriden"
        return
    end
    if !isdir(dirname(fname))
        mkpath(dirname(fname))
    end

    if !testrun
        general_m3e_scan(N, observables, fname, obs_names;
            r1s=vcat(0.0, 10 .^ LinRange(-3, 3, 5)),
            r2s=vcat(0.0, 10 .^ LinRange(-3, 3, 5)),
            r3s=vcat(0.0, 10 .^ LinRange(-3, 3, 5)),
            ets=LinRange(0.0, 6.0, 8),
            ders=LinRange(0.0, 6.0, 8),
            alphas=[0, 0.5],
            ebs=LinRange(0, 6, 3),
            cPs=10 .^ LinRange(-3, 3, 11),
            cRs=10 .^ LinRange(-3, 3, 11),
            do_ss=true
        )
    else
        general_m3e_scan(N, observables, fname, obs_names;
            r1s=10 .^ LinRange(-5, 3, 10),
            alphas=[0, 0.5],
            ets=[0.0, 1.0, 3.0],
            ders=[0.0, 1.0, 3.0],
            ebs=[0.0, 3.0],
            do_ss=true
        )
    end
end

function macroscans(; kwargs...)
    for N in 2:5
        macroscan(N; kwargs...)
    end
end
