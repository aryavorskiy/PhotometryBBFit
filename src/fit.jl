struct PlanckModel end

sparams(::PlanckModel) = (1e10, 5000.0)
min_constraints(::PlanckModel) = (1e7, 1000.0)
max_constraints(::PlanckModel) = (1e70, 1000_000.0)
nparams(model) = length(sparams(model))
spectrum(::PlanckModel, params) = PlanckSpectrum(params...)

using LinearAlgebra

function jacobian(spectrum, filter::Filter)
    return filter_reading(l -> gradient(spectrum, l)[1], filter),
        filter_reading(l -> gradient(spectrum, l)[2], filter)
end
function jacobian!(J, spectrum, pt::SeriesPoint)
    for i in eachindex(pt.filters)
        J[:, i] .= jacobian(spectrum, pt.filters[i])  ./ pt.errs[i]
    end
end

function levenberg_marquardt(model, pt::SeriesPoint; tol = 1e-8, maxiter=10000, verbose=1)
    J = zeros(nparams(model), length(pt.filters))
    β = collect(sparams(model))
    M = zeros(nparams(model), nparams(model))
    rhs = zero(β)
    function step(β, lambda)
        jacobian!(J, spectrum(model, β), pt)
        mul!(M, J, J')
        for i in 1:nparams(model)
            M[i, i] += lambda
        end
        mul!(rhs, J, (pt.mags - filter_reading(spectrum(model, β), pt)) ./ pt.errs)
        verbose > 2 && begin
            @show M
            @show J
            @show rhs
            @show lambda
            println()
        end

        min.(max.(β + M \ rhs, min_constraints(model)), max_constraints(model))
    end

    jacobian!(J, spectrum(model, β), pt)
    λ = maximum(abs, J)^2 / 1e3
    old_chi2 = chi2(spectrum(model, β), pt)
    new_chi2 = chi2(spectrum(model, step(β, λ)), pt)

    while new_chi2 > old_chi2
        λ *= 5
        new_chi2 = chi2(spectrum(model, step(β, λ)), pt)
    end

    for i in 1:maxiter
        newβ = step(β, λ)
        verbose > 1 && (println("iteration #$i:"); @show β)
        new_chi2 = chi2(spectrum(model, newβ), pt)
        if new_chi2 < old_chi2
            λ *= 1.5
        else
            λ /= 5
        end
        old_chi2 = new_chi2
        verbose > 1 && (@show newβ; @show new_chi2; println())

        d = sum(@. abs(β - newβ) / max(β, 1))
        β = newβ
        if d < tol && i > 10
            verbose > 0 && @info "Converged after $i iterations with χ² = $(new_chi2)"
            return spectrum(model, β)
        end
    end
    error("Did not converge!")
end
