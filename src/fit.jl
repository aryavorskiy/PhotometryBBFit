struct BlackBodyModel end

sparams(::BlackBodyModel) = (1e10, 5000.0)
min_constraints(::BlackBodyModel) = (1e7, 1000.0)
max_constraints(::BlackBodyModel) = (1e70, 1000_000.0)
min_constraints(any, i) = min_constraints(any)[i]
max_constraints(any, i) = max_constraints(any)[i]
nparams(model) = length(sparams(model))

spectrum(::BlackBodyModel, params) = PlanckSpectrum(params...)
split_spectrum(spectrum::PlanckSpectrum) = BlackBodyModel(), (spectrum.R, spectrum.T)

using LinearAlgebra

function jacobian(spectrum, filter::Filter)
    return filter_reading(l -> gradient(spectrum, l), filter)
end
function jacobian!(J, spectrum, pt::SeriesPoint)
    for i in eachindex(pt.filters)
        J[:, i] .= jacobian(spectrum, pt.filters[i])  ./ pt.errs[i]
    end
end

struct LMResult{ST}
    spectrum::ST
    covar::Matrix{Float64}
    chi2::Float64
    converged::Bool
    iters::Int
    npoints::Int
end

Base.iterate(r::LMResult) = (r.spectrum, Val(:covar))
Base.iterate(r::LMResult, ::Val{:covar}) = (r.covar, Val(:chi2))
Base.iterate(r::LMResult, ::Val{:chi2}) = (chi2dof(r), Val(:conv))
Base.iterate(r::LMResult, ::Val{:conv}) = (r.converged, Val(:iters))
Base.iterate(r::LMResult, ::Val{:iters}) = (r.iters, Val(:end))
Base.iterate(::LMResult, ::Val{:end}) = nothing

chi2(r::LMResult) = r.chi2
chi2dof(r::LMResult) = r.chi2 / r.npoints

function levenberg_marquardt(model, pt::SeriesPoint; tol = 1e-8, maxiter=1000, verbose=0)
    β = collect(sparams(model))
    newβ = similar(β)

    fr = filter_reading(spectrum(model, β), pt)
    fr_buffer = similar(fr)

    J = zeros(nparams(model), length(pt.filters))
    M = zeros(nparams(model), nparams(model))
    grad = similar(β)

    function step!(fr, β, lambda)
        # assume fr contains readings at β
        jacobian!(J, spectrum(model, β), pt)
        mul!(M, J, J')
        for i in 1:nparams(model)
            M[i, i] += lambda
        end
        mul!(grad, J, (pt.mags .- fr) ./ pt.errs)
        verbose > 1 && begin
            @show M
            @show J
            @show grad
            @show lambda
            println()
        end
        ldiv!(newβ, lu(M), grad)
        for i in 1:nparams(model)
            newβ[i] = min(max(β[i] + newβ[i], min_constraints(model, i)), max_constraints(model, i))
        end
        return newβ
    end
    _chi2(_fr) = sum(((y1, y2, err),) -> ((y1 - y2) / err)^2, zip(pt.mags, _fr, pt.errs))

    jacobian!(J, spectrum(model, β), pt)
    λ₀ = maximum(abs, J)^2 * 1e10
    λ = λ₀
    old_chi2 = _chi2(fr)
    step!(fr, β, λ)
    filter_reading!(fr_buffer, spectrum(model, β), pt)
    new_chi2 = _chi2(fr_buffer)

    converged = false
    iters = 1
    while iters < maxiter && !converged
        # check convergence
        d = sum(@. abs(β - newβ) / max(β, 1))
        if d < tol
            converged = true
            break
        end

        # adjust damping
        verbose > 0 && (println("iteration #$iters:"); @show β)
        rho = (old_chi2 - new_chi2) / sum(abs2, J' * (newβ - β))
        if rho < 5e-2 && λ < λ₀ * 1e16
            # Bad step, increase damping and repeat
            λ *= 10
        else
            if rho > 0.65 && λ > λ₀ * 1e-16
                λ *= 0.2
            end
            old_chi2 = new_chi2
            verbose > 0 && (@show rho; @show newβ; @show new_chi2; println())

            # prepare for new step
            β, newβ = newβ, β
            fr, fr_buffer = fr_buffer, fr
        end
        step!(fr, β, λ)
        filter_reading!(fr_buffer, spectrum(model, newβ), pt)
        new_chi2 = _chi2(fr_buffer)
        iters += 1
    end
    out_spec = spectrum(model, newβ)
    jacobian!(J, out_spec, pt)
    return LMResult(
        out_spec,
        inv(J * J'),
        new_chi2,
        converged,
        iters,
        length(pt.filters))
end
