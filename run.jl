include("src/main.jl")
using Plots
using ProgressMeter
import UnitfulAstro: ly

FILTER_TYPE = :photon

ser = read_series_filterdata("data_13dqy_formatted_for_package.txt", unit=Flux(16e10ly)) do filter
    if isfile("Filters/$filter.txt")
        return read_filter("Filters/$filter.txt", FILTER_TYPE, filter)
    elseif isfile("Filters/$filter.rtf")
        return read_filter("Filters/$filter.rtf", FILTER_TYPE, filter)
    end
end

function plot_RT(series::SeriesFilterdata, model, ts_in=time_domain(series); chi2dof_threshold=100)
    ts = Float64[]
    Rs = Float64[]
    Rerrs = Float64[]
    Ts = Float64[]
    Terrs = Float64[]
    ts_nconverged = Float64[]
    @showprogress for t in ts_in
        res = levenberg_marquardt(model, series(t), verbose=0, tol=1e-8)
        if !res.converged
            push!(ts_nconverged, t)
            continue
        end
        (R, T), covar, χ²dof = res
        (χ²dof > chi2dof_threshold) && continue
        push!(ts, t)
        push!(Rs, R)
        push!(Rerrs, sqrt(covar[1, 1]))
        push!(Ts, T)
        push!(Terrs, sqrt(covar[2, 2]))
    end

    !isempty(ts_nconverged) && @show ts_nconverged
    # pyplot()
    gr()
    p = plot(size=(1000, 450), layout=2, xlab="timestamp, JD", leg=:none)
    scatter!(p[1], ts, Rs, err=Rerrs, ylab = "radius, cm", xlims=(0, NaN), ylims=(0, NaN))
    scatter!(p[2], ts, Ts, err=Terrs, ylab = "temperature, K", xlims=(0, NaN), ylims=(0, 20000))
end

function heatmap_chi2(pt::SeriesPoint, model)
    (R, T), covar, χ²dof = levenberg_marquardt(model, pt)
    d = 0.06
    dT = 400
    local Rs = (1-d)*R:0.02d * R:(1 + d)R
    local Ts = T-dT:dT/50:T+dT
    data = [chi2dof(spectrum(model, (r, t)), pt) for t in Ts, r in Rs]
    # pythonplot()
    gr()
    chistep = √(chi2dof(spectrum(model, (R, T + dT/10)), pt) - χ²dof)
    heatmap(Rs, Ts, data, cbar=true,
        ylab="temperature, K", xlab="radius, cm",
        title="χ²/dof in parameter space", c=:imola)
    contour!(Rs, Ts, data, levels=(χ²dof .+ range(0, step=chistep, length=20) .^ 2), clabels=true, c=:pink)
    scatter!([(R, T)], xerr=[sqrt(covar[1, 1])], yerr=[sqrt(covar[2, 2])], lab="($(trunc(R, sigdigits=3)), $(round(Int, T))), χ²/dof = $(trunc(χ²dof, sigdigits=5))", c=:pink)
end

# heatmap_chi2(ser(10), BlackBodyModel())

dates = [
    9.5866667e-01
    1.1586667e+00
    1.3586667e+00
    1.5586667e+00
    1.7586667e+00
    1.9586667e+00
    2.1586667e+00
    3.1586667e+00
    4.1586667e+00
    5.1586667e+00
    6.1586667e+00
    7.1586667e+00
    8.1586667e+00
    9.1586667e+00
    1.0158667e+01
    1.1158667e+01
    1.6158667e+01
    2.1158667e+01
    2.6158667e+01
    3.1158667e+01
    3.6158667e+01
    4.1158667e+01
    4.6158667e+01
    5.1158667e+01
    5.6158667e+01
]
plot_RT(ser, BlackBodyModel(), dates)
