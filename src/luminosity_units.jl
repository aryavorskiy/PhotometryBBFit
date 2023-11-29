import Unitful: upreferred
import UnitfulAstro: cm
abstract type LightnessUnit end
struct Luminosity <: LightnessUnit end
struct Flux{DistT, AreaT} <: LightnessUnit
    dist::DistT
    area::AreaT
end
Flux(dist) = Flux(dist, cm^2)

"""
    to_luminosity(unit, val, err)
    to_luminosity(unit, series)

converts all units to spectral luminosity in erg/s/AA
"""
to_luminosity(::Luminosity, val, err) = val, err
function to_luminosity(unit::Flux, val, err)
    return (val, err) .* (4pi * upreferred(unit.dist^2 / unit.area))
end
function to_luminosity(unit, ser::Series)
    new_ser = Series(ser.time, similar(ser.mag), similar(ser.err))
    for i in eachindex(ser.time)
        new_ser.mag[i], new_ser.err[i] = to_luminosity(unit, ser.mag[i], ser.err[i])
    end
    return new_ser
end
