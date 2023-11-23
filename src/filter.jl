import Base.Broadcast: broadcasted

struct Filter{T}
    wavelength::Vector{T}
    wavelength_weights::Vector{T}
    transmission::Vector{T}
    norm_const::Float64
    mode::Symbol
    id::String
    function Filter{T}(wavelength::Vector{T}, transmission::Vector{T}, mode::Symbol=:photon, id="") where T
        @assert length(wavelength) == length(transmission)
        @assert mode in (:photon, :energy) "unsupported mode `$mode`"
        wavelength_weights = zero(wavelength)
        wavelength_weights[1:end-1] += diff(wavelength)
        wavelength_weights[2:end] += diff(wavelength)
        wavelength_weights /= 2
        norm_const = mode == :energy ?
        sum(@. wavelength_weights * transmission) :
        sum(@. wavelength_weights * wavelength * transmission)
        return new{T}(wavelength, wavelength_weights, transmission, norm_const, mode, id)
    end
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

using Downloads, Logging

"""
    download_filter(id)

Download filter info from svo2.cab.inta-csic.es website. Returns a `Filter`.
"""
function download_filter(id::String, mode=:photon)
    file = tempname()
    Downloads.download("http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=$id", file)
    filter = read_filter(file, mode, id)
    rm(file)
    return filter
end

"""
    read_filter(file)

Reads filter info from `file` in tab-separated format. Returns a `Filter`.
"""
function read_filter(file, mode=:photon, id=file)
    wavelength = Float64[]
    transparency = Float64[]
    for line in eachline(file)
        wl, ts = parse.(Float64, split(line, r"\s+|,\s*"))
        push!(wavelength, wl)
        push!(transparency, ts)
    end
    length(wavelength) == 0 && @warn "Empty filter data in filename `$file`"
    return Filter{Float64}(wavelength, transparency, mode, id)
end

const c = 3e10
function filter_reading(spectrum, filter::Filter)
    if filter.mode == :energy
        return sum(broadcasted(*, filter.wavelength_weights, broadcasted(spectrum, filter.wavelength), filter.transmission)) / filter.norm_const
    elseif filter.mode == :photon
        return sum(broadcasted(*, filter.wavelength_weights, broadcasted(spectrum, filter.wavelength), filter.transmission, filter.wavelength)) / filter.norm_const
    end
end
