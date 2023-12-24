# cgs units
using StaticArrays
const h = 6.626e-27
const kb = 1.38e-16
const c = 3e10

struct PlanckSpectrum <: AbstractSpectrum
    R::Float64
    T::Float64
end

params(ps::PlanckSpectrum) = (ps.R, ps.T)
params_str(ps::PlanckSpectrum) =
    join(["$k = $(trunc(v, sigdigits=5))" for (k, v) in zip(param_names(ps), params(ps))], ", ")

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

function guess(::BlackBodyModel, pt::SeriesPoint)
    i = findmax(pt.vals)[2]
    T = (h * c * 1e8 / kb) / lambda_eff(pt.filters[i])
    (1e10, T)
end
constraints(::BlackBodyModel) = (1e7 => 1e70, 1e3 => 1e6)
spectrum(::BlackBodyModel, params) = PlanckSpectrum(params...)

abstract type SpectrumWrapper <: AbstractSpectrum end
params(rs::SpectrumWrapper) = params(rs.spectrum)
params_str(rs::SpectrumWrapper) = params_str(rs.spectrum)
param_names(rs::SpectrumWrapper) = param_names(rs.spectrum)
abstract type ModelWrapper <: AbstractModel end
guess(rm::ModelWrapper, pt::SeriesPoint) = guess(rm.model, pt)
constraints(rm::ModelWrapper) = constraints(rm.model)

struct RedshiftSpectrum{ST} <: SpectrumWrapper
    spectrum::ST
    redshift::Float64
end

spectral_density(rs::RedshiftSpectrum, λ) = spectral_density(rs.spectrum, λ / (1 + rs.redshift)) * (1 + rs.redshift)
gradient(rs::RedshiftSpectrum, λ) = gradient(rs.spectrum, λ / (1 + rs.redshift)) * (1 + rs.redshift)

struct RedshiftModel{MT} <: ModelWrapper
    model::MT
    redshift::Float64
end

spectrum(rm::RedshiftModel, params) = RedshiftSpectrum(spectrum(rm.model, params), rm.redshift)

Redshift(model::AbstractModel, z) = RedshiftModel(model, Float64(z))
Redshift(spectrum::AbstractSpectrum, z) = RedshiftSpectrum(spectrum, Float64(z))

struct ExtinctionSpectrum{ST, ET} <: SpectrumWrapper
    spectrum::ST
    extinction::ET
    Av::Float64
end
spectral_density(es::ExtinctionSpectrum, λ) = spectral_density(es.spectrum, λ) * 10^(-0.4 * es.Av * es.extinction(λ))
gradient(es::ExtinctionSpectrum, λ) = gradient(es.spectrum, λ) * 10^(-0.4 * es.Av * es.extinction(λ))

struct ExtinctionModel{MT, ET} <: ModelWrapper
    model::MT
    extinction::ET
    Av::Float64
end

spectrum(em::ExtinctionModel, params) = ExtinctionSpectrum(spectrum(em.model, params), em.extinction, em.Av)

Extinction(model::AbstractModel, ext; Av=1) = ExtinctionModel(model, ext, Float64(Av))
Extinction(spectrum::AbstractSpectrum, ext; Av=1) = ExtinctionSpectrum(spectrum, ext, Float64(Av))
