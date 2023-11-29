# cgs units
using StaticArrays
const h = 6.626e-27
const kb = 1.38e-16
const c = 3e10

struct PlanckSpectrum <: AbstractSpectrum
    R::Float64
    T::Float64
end

Base.iterate(spectrum::PlanckSpectrum) = spectrum.R, Val(:T)
Base.iterate(spectrum::PlanckSpectrum, ::Val{:T}) = spectrum.T, Val(:end)
Base.iterate(::PlanckSpectrum, ::Val{:end}) = nothing
Base.length(::PlanckSpectrum) = 2

@inline spectral_density(spectrum::PlanckSpectrum, λ) =
    (2π * h * c^2 * 1e40 * 1e-8 * 4pi) / λ^5 /
        (exp((h * c * 1e8 / kb) / λ / spectrum.T) - 1) * spectrum.R^2
@inline function gradient(spectrum::PlanckSpectrum, λ)
    pw = (h * c * 1e8 / kb) / λ / spectrum.T
    ex = exp(pw)
    orig = spectral_density(spectrum, λ)
    return SA[orig * 2 / spectrum.R,
    orig / (ex - 1) * ex * pw / spectrum.T]
end

struct BlackBodyModel <: AbstractModel end

start_params(::BlackBodyModel) = (1e10, 5000.0)
constraints(::BlackBodyModel) = (1e7 => 1e70, 1e3 => 1e6)

spectrum(::BlackBodyModel, params) = PlanckSpectrum(params...)
