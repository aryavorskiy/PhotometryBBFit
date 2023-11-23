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

function plot_RT(series::SeriesFilterdata, model)
    ts = Float64[]
    Rs = Float64[]
    Ts = Float64[]
    @showprogress for t in time_domain(series)
        (R, T), chi2 = levenberg_marquardt(model, series(t), verbose=0, tol=1e-10)
        (chi2 > 300) && continue
        push!(ts, t)
        push!(Rs, R)
        push!(Ts, T)
    end

    pyplot()
    p = plot(size=(1000, 450), layout=2, xlab="timestamp, JD", leg=:none)
    scatter!(p[1], ts, Rs, ylab = "radius, cm", xlims=(0, NaN), ylims=(0, NaN))
    scatter!(p[2], ts, Ts, ylab = "temperature, K", xlims=(0, NaN), ylims=(0, 20000))
end

function heatmap_chi2(pt::SeriesPoint, model)
    (R, T), χ² = levenberg_marquardt(model, pt)
    @show χ²
    d = 0.6
    dT = 6000
    local Rs = (1-d)*R:0.02d * R:(1 + d)R
    local Ts = T-dT:dT/50:T+dT
    chi2(spectrum(model, (4, 4)), pt)
    data = [chi2(spectrum(model, (r, t)), pt)
        for t in Ts, r in Rs]
    # pyplot()
    gr()
    chistep = √(chi2(spectrum(model, (R, T + dT/10)), pt) - χ²)
    heatmap(Rs, Ts, data, cbar=true,
        ylab="temperature, K", xlab="radius, cm",
        title="χ² in parameter space", c=:imola)
    contour!(Rs, Ts, data, levels=(χ² .+ range(0, step=chistep, length=20) .^ 2), clabels=true, c=:pink)
    scatter!([(R, T)], lab="($(trunc(R, sigdigits=3)), $(round(Int, T))), χ² = $(trunc(χ², sigdigits=5))", c=:pink)
end

# plot_RT(ser, BlackBodyModel())
heatmap_chi2(ser(10), BlackBodyModel())
