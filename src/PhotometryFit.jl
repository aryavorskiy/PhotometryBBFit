module PhotometryFit

include("core/filter.jl")
export filter_flux, interpolate!, interpolate, lambda_eff
include("core/series.jl")
export fit_summary, time_domain, chi2, chi2dof
include("core/fit.jl")
export fit, chi2dof_threshold, spectrum
include("core/planck.jl")
export BlackBodyModel, PlanckSpectrum, Redshift, Extinction

include("luminosity_units.jl")
export Luminosity, Flux

include("load_utils.jl")
export read_filter, read_photometry_data, FilterFolder, FilterInfo

include("plot_utils.jl")

end # module PhotometryFit
