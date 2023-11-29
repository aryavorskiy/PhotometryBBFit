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

struct SeriesPoint{T}
    filters::Vector{Filter{T}}
    mags::Vector{T}
    errs::Vector{T}
end
function filter_flux(spectrum::AbstractSpectrum, pt::SeriesPoint)
    return [filter_flux(spectrum, f) for f in pt.filters]
end
function filter_flux!(buffer, spectrum::AbstractSpectrum, pt::SeriesPoint)
    for (i, f) in enumerate(pt.filters)
        buffer[i] = filter_flux(spectrum, f)
    end
end
function chi2(spectrum::AbstractSpectrum, pt::SeriesPoint)
    return sum(((y1, y2, err),) -> ((y1 - y2) / err)^2,
        zip(pt.mags, filter_flux(spectrum, pt), pt.errs))
end
chi2dof(spectrum::AbstractSpectrum, pt::SeriesPoint) = chi2(spectrum, pt) / length(pt.filters)
const REPORT_3σ = ("out of 3σ", "in 3σ")
function summary(spectrum::AbstractSpectrum, pt::SeriesPoint)
    println("χ² = $(chi2(spectrum, pt))")
    m = filter_flux(spectrum, pt)
    for i in eachindex(pt.filters)
        println("""$(pt.filters[i].id) (#$i): value $(pt.mags[i]), error $(pt.errs[i])
        Modelled $(m[i]), ratio $(m[i] / pt.mags[i]), $(REPORT_3σ[(abs(pt.mags[i] - m[i]) < 3 * pt.errs[i]) + 1])""")
    end
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
