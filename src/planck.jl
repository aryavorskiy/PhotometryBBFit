# cgs units
using StaticArrays
const h = 6.626e-27
const kb = 1.38e-16

struct PlanckSpectrum
    R::Float64
    T::Float64
end

Base.iterate(spectrum::PlanckSpectrum) = spectrum.R, Val(:T)
Base.iterate(spectrum::PlanckSpectrum, ::Val{:T}) = spectrum.T, Val(:End)
Base.iterate(::PlanckSpectrum, ::Val{:End}) = nothing

@inline (spectrum::PlanckSpectrum)(λ) =
    (2h * c^2 * 1e40 * 1e-8 * 4pi) / λ^5 / (exp((h * c * 1e8 / kb) / λ / spectrum.T) - 1) * spectrum.R^2
@inline function gradient(spectrum::PlanckSpectrum, λ)
    pw = (h * c * 1e8 / kb) / λ / spectrum.T
    ex = exp(pw)
    orig = spectrum(λ)
    return SA[orig * 2 / spectrum.R,
    orig / (ex - 1) * ex * pw / spectrum.T]
end
