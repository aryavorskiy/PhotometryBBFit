module PhotometryFit

include("core/filter.jl")
export filter_flux, interpolate
include("core/series.jl")
export fit_summary, time_domain
include("core/fit.jl")
export fit, chi2dof_threshold
include("core/planck.jl")
export BlackBodyModel, PlanckSpectrum

include("luminosity_units.jl")
export Luminosity, Flux

include("load_utils.jl")
export read_filter, read_photometry_data, SVO2Website, FilterFolder, ConversionRules

include("plot_utils.jl")

end # module PhotometryFit
