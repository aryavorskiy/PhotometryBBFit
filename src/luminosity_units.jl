import Unitful: upreferred, cm
abstract type LightnessUnit end
"""
    to_luminosity(unit, filter, number)
    to_luminosity(unit, series)

Convert all units to spectral luminosity in erg/s/AA
"""
function to_luminosity(unit::LightnessUnit, f::Filter, ser::Series)
    new_ser = Series(ser.time, similar(ser.val), similar(ser.err))
    for i in eachindex(ser.time)
        conv_val = to_luminosity(unit, f, ser.val[i] ± ser.err[i])
        new_ser.val[i] = conv_val.val
        new_ser.err[i] = conv_val.err
    end
    return new_ser
end

struct Luminosity <: LightnessUnit end
to_luminosity(::Luminosity, ::Filter, val::Number) = val

struct Flux{DistT, AreaT} <: LightnessUnit
    dist::DistT
    area::AreaT
end
Flux(;dist, area=cm^2) = Flux(dist, area)
function to_luminosity(unit::Flux, ::Filter, val::Number)
    return val * 4pi * upreferred(unit.dist^2 / unit.area)
end

struct Magnitude{DistT, AreaT} <: LightnessUnit
    dist::DistT
    area::AreaT
end
Magnitude(;dist, area=cm^2) = Magnitude(dist, area)
function zero_flux(system, band)
    zeroflux = 3691
    corr_mag = if system == :ab
        0.0
    elseif system == :vega
        band == :UVW2 ? 1.73 :
        band == :UVM2 ? 1.69 :
        band == :UVW1 ? 1.51 :
        band == :U ? 1.02 :
        band == :B ? -0.13 :
        band == :V ? -0.01 :
        band == :R ? 0.05 :
        band == :I ? -0.31 :
        error("unsupported $band-band in Vega system")
    end
    return zeroflux * 10^(0.4*corr_mag)
end
function to_luminosity(unit::Magnitude, f::Filter, val::Number)
    flux_jy = zero_flux(f.meta.system, f.meta.band) * 10^(-0.4*val)
    flux = flux_jy * 1e23 * (c / lambda_eff(f)^2)
    return flux * 4pi * upreferred(unit.dist^2 / unit.area)
end
