struct Series{T}
    time::Vector{T}
    mag::Vector{T}
    err::Vector{T}
end
Series() = Series{Float64}([], [], [])
function insert_measurement!(ms::Series, ind, val, err)
    i = searchsortedfirst(ms.time, ind)
    insert!(ms.time, i, ind)
    insert!(ms.mag, i, val)
    insert!(ms.err, i, err)
end

function (ser::Series)(time)
    i = findfirst(>(time), ser.time)
    i in (1, nothing) && error("timestamp out of bounds!")
    lw = (time - ser.time[i-1]) / (ser.time[i] - ser.time[i-1])
    return (ser.mag[i-1] * lw + ser.mag[i] * (1 - lw),
    √(ser.err[i-1]^2 * lw + ser.err[i]^2 * (1 - lw)))
end

function read_series(file; dlm=',')
    tags_and_series = Dict{String, Series{Float64}}()
    for line in eachline(file)
        startswith(line, '#') && continue
        time_s, mag_s, err_s, tag = split(line, dlm)
        time, mag, err = parse.(Float64, (time_s, mag_s, err_s))
        tag in keys(tags_and_series) ||
            (tags_and_series[tag] = Series())
        insert_measurement!(tags_and_series[tag], time, mag, err)
    end
    return tags_and_series
end

struct SeriesPoint{T}
    filters::Vector{Filter{T}}
    mags::Vector{T}
    errs::Vector{T}
end
planck_model(sp::SeriesPoint, (R, T)) = [planck_model(filter, (R, T)) for filter in sp.filters]
function resid2(sp::SeriesPoint, (R, T))
    return abs2.(sp.mags .- planck_model(sp, (R, T)))
end
function chi2(sp::SeriesPoint, (R, T))
    re = resid2(sp, (R, T))
    return sum(@. re / sp.errs^2)
end

# Unit conversion

import UnitfulAstro
import UnitfulAstro: cm
abstract type LightnessUnit end
struct Luminosity <: LightnessUnit end
struct Flux{DistT, AreaT} <: LightnessUnit
    dist::DistT
    area::AreaT
end
Flux(dist) = Flux(dist, cm^2)

"""
    convert_unit(val, unit)

converts all units to luminosity
"""
to_luminosity(::Luminosity, val, err) = val, err
function to_luminosity(unit::Flux, val, err)
    return (val, err) .* (4pi * UnitfulAstro.Unitful.upreferred(unit.dist^2 / unit.area))
end
function to_luminosity(unit, ser::Series)
    new_ser = Series(ser.time, similar(ser.mag), similar(ser.err))
    for i in eachindex(ser.time)
        new_ser.mag[i], new_ser.err[i] = to_luminosity(unit, ser.mag[i], ser.err[i])
    end
    return new_ser
end

struct SeriesFilterdata{T}
    filters::Vector{Filter{T}}
    sers::Vector{Series{T}}
end
function (sfd::SeriesFilterdata)(time)
    return SeriesPoint(sfd.filters,
    [ser(time)[1] for ser in sfd.sers],
    [ser(time)[2] for ser in sfd.sers])
end

find_filter(f::Function, tag::String) = f(tag)
function read_series_filterdata(callback, file; unit=Luminosity())
    rser = read_series(file)
    filters = [find_filter(callback, tag) for (tag, ser) in rser]
    sers = [to_luminosity(unit, ser) for (tag, ser) in rser]
    notfound_mask = filters .=== nothing
    if any(notfound_mask)
        @warn """some filter tags were not resolved:
        $(join(collect(keys(rser))[notfound_mask], ", "))"""
        filters = filters[.!notfound_mask]
        sers = sers[.!notfound_mask]
    end
    return SeriesFilterdata(convert(Vector{Filter{Float64}}, filters), sers)
end
