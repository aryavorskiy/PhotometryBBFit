import Base.Broadcast: broadcasted

const FILTER_MODES = (:photon, :energy)
check_mode(mode::Symbol) = @assert mode in FILTER_MODES "unsupported mode '$(mode)'\n" *
"expected one of $FILTER_MODES"

const FILTER_BANDS = (
    :vega => (:UVW2, :UVM2, :UVW1, :U, :B, :V, :R, :I),
    :ab => (:u, :g, :r, :i, :z))
function get_system(band::Symbol)
    for (sys, bands) in FILTER_BANDS
        band in bands && return sys
    end
    return :none
end
_lowercase(sym::Symbol) = Symbol(lowercase(string(sym)))
function check_band(sys::Symbol, band::Symbol)
    if sys == :none
        band == :none && return
        error("Band $band has no known photometric system")
    end
    i = findfirst(a -> first(a) == sys, FILTER_BANDS)
    filter_syss = Tuple(first.(FILTER_BANDS))
    @assert i !== nothing """unknown photometric system '$sys'
    expected one of $filter_syss"""
    bands = FILTER_BANDS[i][2]
    @assert band in bands """unknown band '$band' in system '$sys'
    expected one of $bands"""
end

mutable struct FilterMetadata
    id::String
    family::String
    system::Symbol
    band::Symbol
    mode::Symbol
end
function FilterMetadata(;family="", band=:none, system=get_system(band), mode=:photon, id=family * "_" * string(band))
    check_mode(mode)
    check_band(system, band)
    FilterMetadata(id, lowercase(family), system, band, mode)
end

struct Filter{T}
    wavelength::Vector{T}
    transmission::Vector{T}
    weights::Vector{T}
    meta::FilterMetadata
    function Filter{T}(wavelength::Vector{T}, transmission::Vector{T},
            meta=FilterMetadata(), precalculate_weights=true) where T
        @assert length(wavelength) == length(transmission)
        @assert issorted(wavelength)

        weights = similar(wavelength)
        f = new{T}(wavelength, transmission, weights, meta)
        precalculate_weights && update_weights!(f)
        return f
    end
end
function update_weights!(f::Filter)
    resize!(f.weights, length(f.wavelength))
    fill!(f.weights, 0)
    d = diff(f.wavelength)
    f.weights[1:end-1] += d
    f.weights[2:end] += d
    @. f.weights *= 0.5 * f.transmission
    if f.meta.mode === :photon
        @. f.weights *= f.wavelength
    end
    f.weights ./= sum(f.weights)
    return f
end

Base.:(==)(f1::Filter, f2::Filter) =
    f1.wavelength == f2.wavelength && f1.weights == f2.weights
Base.hash(f::Filter) = hash(f.wavelength, hash(f.weights))
function Base.copy(f::Filter{T}) where T
    f2 = Filter{T}(copy(f.wavelength), copy(f.transmission), f.meta, false)
    copyto!(f2.weights, f.weights)
    return f2
end

Base.summary(io::IO, f::Filter) = summary(io, f.meta)
name(f::Filter) = f.meta.id
function Base.show(io::IO, ::MIME"text/plain", filter::Filter)
    if get(io, :compact, false)
        print(io, name(filter), "@ $(trunc(lambda_eff(filter), digits=2)) AA")
    else
        print("""Filter '$(name(filter))' with `$(filter.meta.mode)`-type detector
        λ_eff = $(trunc(lambda_eff(filter), digits=2)) AA, range $(filter.wavelength[1]) .. $(filter.wavelength[end])""")
    end
end
Base.show(io::IO, filter::Filter) =
    if get(io, :compact, false)
        show(io, MIME"text/plain"(), filter)
    else
        Base.show_default(io, filter)
    end

function Base.write(io::IO, filter::Filter{T}) where T
    write(io, length(filter.transmission)) +
    write(io, filter.wavelength) +
    write(io, filter.transmission)
end

function Base.read(io::IO, ::Type{Filter{T}}) where T
    len = read(io, Int)
    wavelength = Array{T}(undef, len)
    transparency = Array{T}(undef, len)
    read!(io, wavelength)
    read!(io, transparency)
    return Filter{T}(wavelength, transparency)
end

filter_flux(spectrum, filter::Filter) =
    sum(broadcasted(*, filter.weights, broadcasted(spectrum, filter.wavelength)))

lambda_eff(filter::Filter) = filter_flux(l -> l^-2, filter) ^ -0.5
# lambda_eff(filter::Filter) = (filter_flux(l -> l, filter) / filter_flux(l -> l^-1, filter)) ^ -0.5

function interpolate!(filter::Filter{T}, wavelengths) where T
    l = length(wavelengths)
    transmissions = T[]
    for wl in wavelengths
        i = findfirst(>(wl), filter.wavelength)
        if i in (1, nothing)
            push!(transmissions, 0)
        else
            lw = (wl - filter.wavelength[i-1]) / (filter.wavelength[i] - filter.wavelength[i-1])
            tr = filter.transmission[i-1] * (1 - lw) + filter.transmission[i] * lw
            push!(transmissions, tr)
        end
    end
    resize!(filter.transmission, l)
    copyto!(filter.transmission, transmissions)
    resize!(filter.wavelength, l)
    copyto!(filter.wavelength, wavelengths)
    update_weights!(filter)
end
interpolate!(filter::Filter; step) =
    interpolate!(filter, minimum(filter.wavelength):step:maximum(filter.wavelength))

function interpolate(filter::Filter{T}, wavelengths) where T
    interpolate!(copy(filter), wavelengths)
end
interpolate(filter::Filter; step) =
    interpolate(filter, minimum(filter.wavelength):step:maximum(filter.wavelength))
