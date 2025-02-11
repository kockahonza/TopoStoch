using DrWatson
@quickactivate "cluster_env"

include("../../scripts/m3_scanner.jl")

function bigrun2(N=4; testrun=false)
    observables = [ObsGoingRoundTheLoop{Matrix,false}(N), ObsGoingRoundTheLoop{Mean,false}(N)]
    obs_names = ["rtl_mats", "rtl_means"]

    fname = f"./bigrun2_N{N}.jld2"
    if ispath(fname)
        @error f"Path {fname} already exists, exiting so that it is not overriden"
    end
    if !isdir(dirname(fname))
        mkpath(dirname(fname))
    end

    if !testrun
        general_m3_scan(N, observables, fname, obs_names;
            r1s=10 .^ LinRange(-5, 3, 15),
            r2s=10 .^ LinRange(-5, 3, 15),
            r3s=10 .^ LinRange(-5, 3, 15),
            alphas=[0, 0.5],
            ebs=LinRange(0, 6, 5),
            cPs=10 .^ LinRange(-5, 1, 15),
            cRs=10 .^ LinRange(-1, 2, 15),
            do_ss=false
        )
    else
        general_m3_scan(N, observables, fname, obs_names;
            r1s=10 .^ LinRange(-5, 3, 15),
            alphas=[0, 0.5],
            ebs=LinRange(0, 6, 5),
            do_ss=false
        )
    end
end

function bigrun2ss(N=4; testrun=false, do_ss=true)
    observables = [ObsGoingRoundTheLoop{Matrix,do_ss}(N), ObsGoingRoundTheLoop{Mean,do_ss}(N)]
    obs_names = ["rtlss_mats", "rtlss_means"]

    fname = f"./bigrun2ss_N{N}.jld2"
    if ispath(fname)
        @error f"Path {fname} already exists, exiting so that it is not overriden"
    end
    if !isdir(dirname(fname))
        mkpath(dirname(fname))
    end

    if !testrun
        general_m3_scan(N, observables, fname, obs_names;
            r1s=10 .^ LinRange(-5, 3, 15),
            r2s=10 .^ LinRange(-5, 3, 15),
            r3s=10 .^ LinRange(-5, 3, 15),
            alphas=[0, 0.5],
            ebs=LinRange(0, 6, 5),
            cPs=10 .^ LinRange(-5, 1, 15),
            cRs=10 .^ LinRange(-1, 2, 15),
            do_ss=do_ss
        )
    else
        general_m3_scan(N, observables, fname, obs_names;
            r1s=10 .^ LinRange(-5, 3, 15),
            alphas=[0, 0.5],
            ebs=LinRange(0, 6, 5),
            do_ss=do_ss
        )
    end
end
