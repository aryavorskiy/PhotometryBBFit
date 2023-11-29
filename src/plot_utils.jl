using RecipesBase

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
        ylims := (0, 20000)
        yformatter := :plain
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

@recipe function f(res::LMResult, ::Val{:heatmap})
    (R, T), covar, χ²dof = res
    model = res.model
    d = 0.06
    dT = 400
    local Rs = (1-d)*R:0.02d * R:(1 + d)R
    local Ts = T-dT:dT/50:T+dT
    data = [chi2dof(spectrum(model, (r, t)), res.pt) for t in Ts, r in Rs]
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
        seriescolor := :darkgrey
        linestyle := :dash
        label := :none
        seriestype := :vline
        (1 - 0.9d)*R:0.3d * R:(1 + 0.9d)R
    end
    @series begin
        seriescolor := :darkgrey
        linestyle := :dash
        label := :none
        seriestype := :hline
        T-0.9dT:dT*0.3:T+0.9dT
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

@recipe function f(res::LMResult)
    seriestype --> :heatmap # TODO change to :path
    if plotattributes[:seriestype] == :heatmap
        (res, Val(:heatmap))
    elseif plotattributes[:seriestype] == :path
        (res, Val(:sed))    # TODO
    else error("unsupported series type $(plotattributes[:seriestype])")
    end
end
