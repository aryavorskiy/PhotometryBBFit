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
    checkbounds(Bool, ser, time) || error("timestamp out of bounds!")
    i = findfirst(≥(time), ser.time)
    i === nothing && error("timestamp $time out of bounds $(extrema(ser.time))!")
    i == 1 && return ser.mag[1], ser.err[1]
    lw = (time - ser.time[i-1]) / (ser.time[i] - ser.time[i-1])
    return (ser.mag[i-1] * lw + ser.mag[i] * (1 - lw),
    √(ser.err[i-1]^2 * lw + ser.err[i]^2 * (1 - lw)))
end

Base.checkbounds(::Type{Bool}, ser::Series, time) = ser.time[begin] ≤ time ≤ ser.time[end]

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
    tags = collect(keys(tags_and_series))
    sers = collect(values(tags_and_series))
    sp = sortperm(tags)
    return tags[sp], sers[sp]
end

struct SeriesPoint{T}
    filters::Vector{Filter{T}}
    mags::Vector{T}
    errs::Vector{T}
end
function filter_reading(spectrum, pt::SeriesPoint)
    return [filter_reading(spectrum, f) for f in pt.filters]
end
function filter_reading!(buffer, spectrum, pt::SeriesPoint)
    for (i, f) in enumerate(pt.filters)
        buffer[i] = filter_reading(spectrum, f)
    end
end
function chi2(spectrum, pt::SeriesPoint)
    return sum(((y1, y2, err),) -> ((y1 - y2) / err)^2,
        zip(pt.mags, filter_reading(spectrum, pt), pt.errs))
end
const REPORT_3σ = ("out of 3σ", "in 3σ")
function summary(spectrum, pt::SeriesPoint)
    println("χ² = $(chi2(spectrum, pt))")
    m = filter_reading(spectrum, pt)
    for i in eachindex(pt.filters)
        println("""$(pt.filters[i].id) (#$i): value $(pt.mags[i]), error $(pt.errs[i])
        Modelled $(m[i]), ratio $(m[i] / pt.mags[i]), $(REPORT_3σ[(abs(pt.mags[i] - m[i]) < 3 * pt.errs[i]) + 1])""")
    end
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
    mask = checkbounds.(Bool, sfd.sers, time)
    return SeriesPoint(sfd.filters[mask],
    [ser(time)[1] for ser in sfd.sers[mask]],
    [ser(time)[2] for ser in sfd.sers[mask]])
end
function time_domain(sfd::SeriesFilterdata; tol=1e-8, n_overlap=2)
    ts = Float64[]
    tstart = sort([ser.time[1] for ser in sfd.sers])[n_overlap]
    tend = sort([ser.time[end] for ser in sfd.sers])[end - n_overlap + 1]
    for ser in sfd.sers
        for t in ser.time
            tstart ≤ t ≤ tend || continue
            any(tt -> isapprox(t, tt, atol=tol), ts) || push!(ts, t)
        end
    end
    return sort(ts)
end

find_filter(f::Function, tag::String) = f(tag)
function read_series_filterdata(callback, file; unit=Luminosity())
    tags, n_sers = read_series(file)
    filters = [find_filter(callback, tag) for tag in tags]
    sers = [to_luminosity(unit, ser) for ser in n_sers]
    notfound_mask = filters .=== nothing
    if any(notfound_mask)
        @warn """some filter tags were not resolved:
        $(join(collect(keys(rser))[notfound_mask], ", "))"""
        filters = filters[.!notfound_mask]
        sers = sers[.!notfound_mask]
    end
    return SeriesFilterdata(convert(Vector{Filter{Float64}}, filters), sers)
end
