using ProgressMeter
using RecipesBase

struct FitSeries{ST, MT, T}
    timestamps::Vector{T}
    fitresults::Vector{LMResult{ST, MT, T}}
end
function fit(ser::SeriesFilterdata, model, times=time_domain(ser);
        filter=true, threshold=Inf, kw...)
    fser = FitSeries(times,
        @showprogress[levenberg_marquardt(model, ser(t); kw...) for t in times])
    if filter
        return filter_by_chi2dof!(fser; threshold=threshold)
    else return fser
    end
end
function filter_by_chi2dof!(fser::FitSeries; threshold=Inf)
    inds = findall(fser.fitresults) do res
        !res.converged || chi2dof(res) > threshold
    end
    deleteat!(fser.timestamps, inds)
    deleteat!(fser.fitresults, inds)
    return fser
end
filter_by_chi2dof(fser; kw...) =
    filter_by_chi2dof!(FitSeries(copy(fser.timestamps), copy(fser.fitresults)); kw...)

function Base.getindex(fser::FitSeries; t)
    i = findmin(_t -> abs(_t - t), fser.timestamps)[2]
    return fser[i]
end
Base.getindex(fser::FitSeries, i::Int) = fser.fitresults[i]

struct ParamSeries{T}
    xs::Vector{T}
    ys::Vector{T}
    yerrs::Vector{T}
end
@recipe function f(ps::ParamSeries)
    yerror := ps.yerrs
    lims := :round
    ps.xs, ps.ys
end
@recipe function f(fser::FitSeries)
    layout := 2
    xlabel := "timestamp"
    legend := :none
    @series begin
        subplot := 1
        ylabel := "radius, cm"
        fser.R
    end
    @series begin
        subplot := 2
        ylabel := "temperature, K"
        fser.T
    end
end

function Base.getproperty(fser::FitSeries{ST}, param::Symbol) where ST
    param_index = findfirst(==(param), fieldnames(ST))
    if param_index !== nothing
        return ParamSeries(fser.timestamps,
        [getproperty(res.spectrum, param) for res in fser.fitresults],
        [sqrt(res.covar[param_index, param_index]) for res in fser.fitresults])
    end
    return getfield(fser, param)
end

@recipe function f(res::LMResult)
    (R, T), covar, χ²dof = res
    model = res.model
    d = 0.06
    dT = 400
    local Rs = (1-d)*R:0.02d * R:(1 + d)R
    local Ts = T-dT:dT/50:T+dT
    data = [chi2dof(spectrum(model, (r, t)), res.pt) for t in Ts, r in Rs]
    gr()
    @series begin
        seriestype := :heatmap
        title --> "χ²/dof in parameter space"
        seriescolor --> :imola
        (Rs, Ts, data)
    end
    @series begin
        seriestype := :contour
        levels := [chi2dof(spectrum(model, (R * (1 + d * i), T)), res.pt) for i in -1:0.1:1]
        seriescolor := get(plotattributes, :linecolor, :pink)
        (Rs, Ts, data)
    end
    @series begin
        seriestype := :scatter
        xerr := [sqrt(covar[1, 1])]
        yerr := [sqrt(covar[2, 2])]
        label := "($(trunc(R, sigdigits=3)), $(round(Int, T))), χ²/dof = $(trunc(χ²dof, sigdigits=5))"
        seriescolor := get(plotattributes, :linecolor, :pink)
        xlims := ((1-d) * R, (1+d) * R)
        ylims := (T - dT, T + dT)
        [(R, T)]
    end
end
