struct Series{T}
    time::Vector{T}
    val::Vector{T}
    err::Vector{T}
end
Series() = Series{Float64}([], [], [])
function insert_measurement!(ms::Series, ind, val, err)
    i = searchsortedfirst(ms.time, ind)
    insert!(ms.time, i, ind)
    insert!(ms.val, i, val)
    insert!(ms.err, i, err)
end
function Base.merge!(s1::Series, s2::Series)
    append!(s1.time, s2.time)
    append!(s1.val, s2.val)
    append!(s1.err, s2.err)
    p = sortperm(s1.time)
    permute!(s1.time, p)
    permute!(s1.val, p)
    permute!(s1.err, p)
end

function interpolate(ser::Series, time)
    checkbounds(Bool, ser, time) || error("timestamp out of bounds!")
    i = findfirst(≥(time), ser.time)
    i === nothing && error("timestamp $time out of bounds $(extrema(ser.time))!")
    i == 1 && return ser.val[1], ser.err[1]
    lw = (time - ser.time[i-1]) / (ser.time[i] - ser.time[i-1])
    return (ser.val[i-1] * (1 - lw) + ser.val[i] * lw,
    √(ser.err[i-1]^2 * (1 - lw) + ser.err[i]^2 * lw),
    min(time - ser.time[i-1], ser.time[i] - time))
end

Base.checkbounds(::Type{Bool}, ser::Series, time) = ser.time[begin] ≤ time ≤ ser.time[end]

struct SeriesPoint{T}
    filters::Vector{Filter{T}}
    vals::Vector{T}
    errs::Vector{T}
end
function filter_flux(spectrum, pt::SeriesPoint)
    return [filter_flux(spectrum, f) for f in pt.filters]
end
function filter_flux!(buffer, spectrum, pt::SeriesPoint)
    @inbounds @simd for i in eachindex(pt.filters)
        buffer[i] = filter_flux(spectrum, pt.filters[i])
    end
end
function chi2(spectrum, pt::SeriesPoint)
    return sum(((y1, y2, err),) -> ((y1 - y2) / err)^2,
        zip(pt.vals, filter_flux(spectrum, pt), pt.errs))
end
chi2dof(spectrum, pt::SeriesPoint) = chi2(spectrum, pt) / length(pt.filters)
const REPORT_3σ = ("out of 3σ", "in 3σ")
function fit_summary(spectrum, pt::SeriesPoint)
    println("χ² = $(chi2(spectrum, pt))")
    println("χ²/dof = $(chi2dof(spectrum, pt))\n")
    m = filter_flux(spectrum, pt)
    for i in eachindex(pt.filters)
        println("""$(pt.filters[i].id) (#$i): $(trunc(pt.vals[i] ± pt.errs[i], sigdigits=4))
        Modelled $(m[i]), normres $(abs(m[i] - pt.vals[i]) / pt.errs[i])""")
    end
end

struct PhotometryData{T} <: AbstractDict{Filter{T}, Series{T}}
    filters::Vector{Filter{T}}
    sers::Vector{Series{T}}
    function PhotometryData{T}(filters::Vector{Filter{T}}, sers::Vector{Series{T}}) where T
        @assert length(filters) == length(sers)
        @inline function filter_lt(f1, f2)
            t1 = f1.id[1:3]
            t2 = f2.id[1:3]
            return t1 == t2 ? isless(lambda_eff(f1), lambda_eff(f2)) : isless(t1, t2)
        end
        sp = sortperm(filters, lt=filter_lt)
        return new{T}(filters[sp], sers[sp])
    end
end
function Base.get(sfd::PhotometryData, f::Filter, default)
    i = findfirst(==(f), sfd.filters)
    i === nothing && return default
    return sfd.sers[i]
end
function Base.get(sfd::PhotometryData, id::String, default)
    i = findfirst(f -> f.id == id, sfd.filters)
    i === nothing && return default
    return sfd.sers[i]
end
Base.iterate(sfd::PhotometryData, i = 1) =
    i ≤ length(sfd) ? (sfd.filters[i] => sfd.sers[i], i + 1) : nothing
Base.length(sfd::PhotometryData) = length(sfd.filters)

function (sfd::PhotometryData)(time; threshold=Inf)
    mask = [checkbounds(Bool, ser, time) && interpolate(ser, time)[3] < threshold
        for ser in sfd.sers]
    return SeriesPoint(sfd.filters[mask],
    [interpolate(ser, time)[1] for ser in sfd.sers[mask]],
    [interpolate(ser, time)[2] for ser in sfd.sers[mask]])
end
function time_domain(sfd::PhotometryData; tol=1e-8, n_overlap=2)
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
