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

@recipe function f(ser::Series)
    xlabel --> "timestamp"
    yerror := ser.err
    ser.time, ser.val
end

@recipe function f(fser::FitSeries)
    size --> (800, 350)
    layout := 2
    xlabel --> "timestamp"
    legend --> :none
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
    Re = sqrt(covar[1, 1])
    Te = sqrt(covar[2, 2])
    dR = 3Re
    dT = 3Te
    local Rs = R - dR:0.02dR:R + dR
    local Ts = T - dT:0.02dT:T + dT
    data = [chi2dof(spectrum(model, (r, t)), res.pt) for t in Ts, r in Rs]
    @series begin
        seriestype := :heatmap
        title --> "χ²/dof in parameter space"
        seriescolor --> :imola
        (Rs, Ts, data)
    end
    @series begin
        seriestype := :contour
        levels := [chi2dof(spectrum(model, (R + dR * i, T)), res.pt) for i in -1:0.1:1]
        seriescolor := get(plotattributes, :linecolor, :pink)
        (Rs, Ts, data)
    end
    @series begin
        seriescolor := :darkgrey
        linestyle := :dash
        label := :none
        seriestype := :vline
        R-0.9dR:0.3dR:R+0.9dR
    end
    @series begin
        seriescolor := :darkgrey
        linestyle := :dash
        label := :none
        seriestype := :hline
        T-0.9dT:0.3dT:T+0.9dT
    end
    @series begin
        seriestype := :scatter
        xerr := [Re]
        yerr := [Te]
        label := "($(trunc(R, sigdigits=3)), $(round(Int, T))), χ²/dof = $(trunc(χ²dof, sigdigits=5))"
        seriescolor := get(plotattributes, :linecolor, :pink)
        xlims := (R - dR, R + dR)
        ylims := (T - dT, T + dT)
        [(R, T)]
    end
end

@recipe function f(res::LMResult, ::Val{:sed})
    lmin = Inf
    lmax = -Inf
    λeffs = Float64[]
    markers = [:diamond, :ltriangle, :rtriangle, :star5, :star8, :rect]
    marker_i = 0
    prefix = ""
    for (i, filter) in enumerate(res.pt.filters)
        if prefix != filter.id[1:3]
            prefix = filter.id[1:3]
            marker_i += 1
        end
        λeff = lambda_eff(filter)
        push!(λeffs, λeff)
        lmin = min(lmin, minimum(filter.wavelength))
        lmax = max(lmax, maximum(filter.wavelength))
        @series begin
            label := filter.id
            yerror := res.pt.errs[i]
            seriestype := :scatter
            markershape := markers[(marker_i - 1) % length(markers) + 1]
            [(λeff, res.pt.vals[i])]
        end
    end
    (plmin, plmax) = get(plotattributes, :xlims, (lmin, lmax))
    isfinite(plmin) && (lmin = plmin)
    isfinite(plmax) && (lmax = plmax)
    @series begin
        title --> params_str(res.spectrum)
        label --> "fit (χ²/dof = $(trunc(chi2dof(res), digits=2)))"
        xformatter := :plain
        xlims := (lmin, lmax)
        l -> res.spectrum(l)
    end
    get(plotattributes, :showscatter, false) && @series begin
        label := ""
        markershape := :cross
        seriestype := :scatter
        (lambda_eff.(res.pt.filters), filter_flux(res.spectrum, res.pt))
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
