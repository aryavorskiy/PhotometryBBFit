using RecipesBase

@recipe function f(filter::Filter)
    label --> filter.id
    @series begin
        fill := true
        fillalpha := 0.3
        fillstyle := :/
        ylims := (0, NaN)
        (filter.wavelength), filter.transmission
    end
    @series begin
        label := :none
        seriestype := :vline
        linestyle := :dash
        seriescolor := get(plotattributes, :seriescolor, :black)
        [lambda_eff(filter)]
    end
end

@recipe function f(fser::FitSeries)
    layout := 2
    xlabel := "timestamp"
    legend := :none
    @series begin
        subplot := 1
        ylabel := "radius, cm"
        fser.timestamps, fser.R
    end
    @series begin
        subplot := 2
        ylabel := "temperature, K"
        ylims := (0, 20000)
        yformatter := :plain
        fser.timestamps, fser.T
    end
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

@recipe function f(res::LMResult, ::Val{:sed}; showscatter=false)
    lmin = Inf
    lmax = -Inf
    λeffs = Float64[]
    markers = [:diamond, :circle, :ltriangle, :rtriangle, :star5, :star8]
    for (i, filter) in enumerate(res.pt.filters)
        λeff = lambda_eff(filter)
        push!(λeffs, λeff)
        lmin = min(lmin, minimum(filter.wavelength))
        lmax = max(lmax, maximum(filter.wavelength))
        @series begin
            label := filter.id
            yerror := res.pt.errs[i]
            seriestype := :scatter
            markershape := markers[(i - 1) % length(markers) + 1]
            [(λeff, res.pt.vals[i])]
        end
    end
    @series begin
        label := ""
        range(lmin, lmax, 100), l -> res.spectrum(l)
    end
    if showscatter
        @series begin
            label := ""
            seriestype := :scatter
            λeffs, [filter_flux(res.spectrum, f) for f in res.pt.filters]
        end
    end
end

@recipe function f(res::LMResult)
    seriestype --> :path
    if plotattributes[:seriestype] == :heatmap
        (res, Val(:heatmap))
    elseif plotattributes[:seriestype] == :path
        (res, Val(:sed))
    else error("unsupported series type $(plotattributes[:seriestype])")
    end
end
