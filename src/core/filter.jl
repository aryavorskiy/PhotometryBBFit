import Base.Broadcast: broadcasted

struct Filter{T}
    wavelength::Vector{T}
    wavelength_weights::Vector{T}
    transmission::Vector{T}
    mode::Symbol
    id::String
    function Filter{T}(wavelength::Vector{T}, transmission::Vector{T}, mode::Symbol=:photon, id="") where T
        @assert length(wavelength) == length(transmission)
        @assert issorted(wavelength)
        @assert mode in (:photon, :energy) "unsupported mode `$mode`"
        wavelength_weights = zero(wavelength)
        wavelength_weights[1:end-1] += diff(wavelength)
        wavelength_weights[2:end] += diff(wavelength)
        wavelength_weights /= 2
        norm_const = mode == :energy ?
        sum(@. wavelength_weights * transmission) :
        sum(@. wavelength_weights * wavelength * transmission)
        transmission ./= norm_const
        return new{T}(wavelength, wavelength_weights, transmission, mode, id)
    end
end

Base.:(==)(f1::Filter, f2::Filter) =
    f1.wavelength == f2.wavelength && f1.transmission == f2.transmission && f1.mode == f2.mode
Base.hash(f::Filter) =
    Base.hash(f.wavelength, hash(f.transmission, hash(f.mode)))

function Base.show(io::IO, ::MIME"text/plain", filter::Filter)
    print(io, "Filter `$(filter.id)`")
    if !get(io, :compact, false)
        print(""" with `$(filter.mode)`-type detector
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

function filter_flux(spectrum, filter::Filter)
    if filter.mode == :energy
        return sum(broadcasted(*, filter.wavelength_weights, broadcasted(spectrum, filter.wavelength), filter.transmission))
    elseif filter.mode == :photon
        return sum(broadcasted(*, filter.wavelength_weights, broadcasted(spectrum, filter.wavelength), filter.transmission, filter.wavelength))
    end
end

lambda_eff(filter::Filter) = filter_flux(l -> l^-2, filter) ^ -0.5

function interpolate(filter::Filter{T}, wavelengths) where T
    transmissions = T[]
    for wl in wavelengths
        i = findfirst(>(wl), filter.wavelength)
        if i in (1, nothing)
            push!(transmissions, 0)
        else
            lw = (wl - filter.wavelength[i-1]) / (filter.wavelength[i] - filter.wavelength[i-1])
            tr = filter.wavelength[i-1] * lw + filter.wavelength[i] * (1 - lw)
            push!(transmissions, tr)
        end
    end
    Filter{T}(wavelengths, transmissions, filter.mode, filter.id)
end
interpolate(filter::Filter; step) =
    interpolate(filter, minimum(filter.wavelength):step:maximum(filter.wavelength))
