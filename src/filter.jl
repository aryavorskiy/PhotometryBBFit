struct Filter{T}
    wavelength::Vector{T}
    wavelength_weights::Vector{T}
    transparency::Vector{T}
    function Filter{T}(wavelength::Vector{T}, transparency::Vector{T}) where T
        @assert length(wavelength) == length(transparency)
        wavelength_weights = zero(wavelength)
        wavelength_weights[1:end-1] += diff(wavelength)
        wavelength_weights[2:end] += diff(wavelength)
        wavelength_weights ./= 2
        return new{T}(wavelength, wavelength_weights, transparency)
    end
end

function Base.write(io::IO, filter::Filter{T}) where T
    write(io, length(filter.transparency)) +
    write(io, filter.wavelength) +
    write(io, filter.transparency)
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
function download_filter(id::String)
    file = tempname()
    Downloads.download("http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=$id", file)
    filter = read_filter(file)
    rm(file)
    return filter
end

"""
    read_filter(file)

Reads filter info from `file` in tab-separated format. Returns a `Filter`.
"""
function read_filter(file)
    wavelength = Float64[]
    transparency = Float64[]
    for line in eachline(file)
        wl, ts = parse.(Float64, split(line, r"\s+|,\s*"))
        push!(wavelength, wl)
        push!(transparency, ts)
    end
    length(wavelength) == 0 && @warn "Empty filter data in filename `$file`"
    return Filter{Float64}(wavelength, transparency)
end

# cgs units
const c = 3e10
const h = 6.626e-27
const kb = 1.38e-16

planck_f(λ, T) = (2h * c^2 * 1e40) / λ^5 / (exp((h * c * 1e8 / kb) / λ / T) - 1)
function planck_model(filter::Filter, (R, T))
    return sum(@. filter.wavelength_weights * planck_f(filter.wavelength, T) * filter.transparency) * 4pi * R^2
end
