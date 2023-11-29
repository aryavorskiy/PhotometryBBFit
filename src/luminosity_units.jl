import Unitful: upreferred, cm
abstract type LightnessUnit end
struct Luminosity <: LightnessUnit end
struct Flux{DistT, AreaT} <: LightnessUnit
    dist::DistT
    area::AreaT
end
Flux(;dist, area=cm^2) = Flux(dist, area)

"""
    to_luminosity(unit, val, err)
    to_luminosity(unit, series)

Convert all units to spectral luminosity in erg/s/AA
"""
to_luminosity(::Luminosity, val::Number) = val
function to_luminosity(unit::Flux, val::Number)
    return val * 4pi * upreferred(unit.dist^2 / unit.area)
end
function to_luminosity(unit::LightnessUnit, ser::Series)
    new_ser = Series(ser.time, similar(ser.val), similar(ser.err))
    for i in eachindex(ser.time)
        conv_val = to_luminosity(unit, ser.val[i] ± ser.err[i])
        new_ser.val[i] = conv_val.val
        new_ser.err[i] = conv_val.err
    end
    return new_ser
end
