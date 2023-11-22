# cgs units
const h = 6.626e-27
const kb = 1.38e-16

struct PlanckModel
    R::Float64
    T::Float64
end


(pm::PlanckModel)(λ) =
    (2h * c^2 * 1e40) / λ^5 / (exp((h * c * 1e8 / kb) / λ / pm.T) - 1) * 1e-8 * 4pi * pm.R^2
function planck_model(filter::Filter, (R, T))

end

planck_df(λ, T) = -(2h^2 * c^3 * 1e48) / λ^6 / (exp((h * c * 1e8 / kb) / λ / T) - 1)^2 * exp((h * c * 1e8 / kb) / λ / T) / kb / T^2
function planck_model(filter::Filter, (R, T), ::Val{:R})
    return sum(@. filter.wavelength_weights * planck_f(filter.wavelength, T) * filter.transmission) * 8pi * R
end
function planck_model(filter::Filter, (R, T), ::Val{:T})
    return sum(@. filter.wavelength_weights * planck_df(filter.wavelength, T) * filter.transmission) * 4pi * R^2
end

planck_model(sp::SeriesPoint, (R, T)) = [planck_model(filter, (R, T)) for filter in sp.filters]
