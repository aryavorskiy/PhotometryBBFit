struct BlackBodyModel end

sparams(::BlackBodyModel) = (1e10, 5000.0)
min_constraints(::BlackBodyModel) = (1e7, 1000.0)
max_constraints(::BlackBodyModel) = (1e70, 1000_000.0)
min_constraints(any, i) = min_constraints(any)[i]
max_constraints(any, i) = max_constraints(any)[i]
nparams(model) = length(sparams(model))

spectrum(::BlackBodyModel, params) = PlanckSpectrum(params...)
split_spectrum(spectrum::PlanckSpectrum) = PlanckModel(), (spectrum.R, spectrum.T)

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
    newβ = zero(β)
    function step!(β, lambda)
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
        ldiv!(newβ, lu(M), rhs)
        for i in 1:nparams(model)
            newβ[i] = min(max(β[i] + newβ[i], min_constraints(model, i)), max_constraints(model, i))
        end
        return newβ
    end

    jacobian!(J, spectrum(model, β), pt)
    λ = maximum(abs, J)^2 / 1e3
    old_chi2 = chi2(spectrum(model, β), pt)
    new_chi2 = chi2(spectrum(model, step!(β, λ)), pt)

    # while new_chi2 > old_chi2
    #     λ *= 5
    #     new_chi2 = chi2(spectrum(model, step!(β, λ)), pt)
    # end

    for i in 1:maxiter
        step!(β, λ)
        verbose > 1 && (println("iteration #$i:"); @show β)
        new_chi2 = chi2(spectrum(model, newβ), pt)
        if new_chi2 < old_chi2
            λ *= 1.5
        elseif λ > eps()
            λ /= 5
        end
        old_chi2 = new_chi2
        verbose > 1 && (@show newβ; @show new_chi2; println())

        d = sum(@. abs(β - newβ) / max(β, 1))
        β, newβ = newβ, β
        if d < tol && i > 10
            verbose > 0 && @info "Converged after $i iterations with χ² = $(new_chi2)"
            return spectrum(model, β), new_chi2
        end
    end
    error("Did not converge!")
end
# function findmin_newt(f; maxiter=10000, tol=1e-2, arg=100., step=1, verbose=false)
#     v = f(arg)
#     for i in 1:maxiter
#         oarg = arg
#         f1 = f(arg + step)
#         f2 = f(arg + 2step)
#         darg = (4 * f1 - f2 - 3 * v) / 3 / abs(f2 - 2 * f1 + v) * step
#         if isfinite(darg)
#             arg -= darg
#             v = f(arg)
#         end
#         verbose && i % 100 == 0 && @show (arg, abs(oarg-arg), v)

#         if abs(oarg - arg) / (arg + abs(oarg - arg)) < tol
#             verbose && @info "Converged on iteration #$i"
#             return arg
#         end
#     end
#     error("Did not converge!")
# end

# using AffineInvariantMCMC

# function find_RT(sp::SeriesPoint; nwalkers=100)
#     numdims = 5
#     thinning = 10
#     numsamples_perwalker = 1000
#     burnin = 100

#     stds = exp(5 * randn(numdims))
#     means = 1 + 5 * rand(numdims)
#     llhood = x->begin
#         retval = 0.
#         for i = eachindex(x)
#             retval -= .5 * ((x[i] - means[i]) / stds[i]) ^ 2
#         end
#         return retval
#     end
#     x0 = rand(numdims, numwalkers) * 10 - 5
#     chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, x0, burnin, 1)
#     chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, chain[:, :, end], numsamples_perwalker, thinning)
#     flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)
# end
