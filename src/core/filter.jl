import Base.Broadcast: broadcasted

struct Filter{T}
    wavelength::Vector{T}
    transmission::Vector{T}
    weights::Vector{T}
    mode::Symbol
    id::String
    function Filter{T}(wavelength::Vector{T}, transmission::Vector{T},
            mode::Symbol=:photon, id="", precalculate_weights=true) where T
        @assert length(wavelength) == length(transmission)
        @assert issorted(wavelength)
        @assert mode in (:photon, :energy) "unsupported mode `$mode`"
        weights = similar(wavelength)
        f = new{T}(wavelength, transmission, weights, mode, id)
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
    if f.mode === :photon
        @. f.weights *= f.wavelength
    end
    f.weights ./= sum(f.weights)
    return f
end

Base.:(==)(f1::Filter, f2::Filter) =
    f1.wavelength == f2.wavelength && f1.transmission == f2.transmission && f1.mode == f2.mode
Base.hash(f::Filter) =
    Base.hash(f.wavelength, hash(f.transmission, hash(f.mode)))
function Base.copy(f::Filter{T}) where T
    f2 = Filter{T}(copy(f.wavelength), copy(f.transmission), f.mode, f.id, false)
    copyto!(f2.weights, f.weights)
    return f2
end

function Base.show(io::IO, ::MIME"text/plain", filter::Filter)
    if get(io, :compact, false)
        print(io, "`$(filter.id)` @ $(trunc(lambda_eff(filter), digits=2)) AA")
    else
        print("""Filter $(filter.id) with `$(filter.mode)`-type detector
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

function interpolate!(filter::Filter{T}, wavelengths) where T
    l = length(wavelengths)
    transmissions = T[]
    for wl in wavelengths
        i = findfirst(>(wl), filter.wavelength)
        if i in (1, nothing)
            push!(transmissions, 0)
        else
            lw = (wl - filter.wavelength[i-1]) / (filter.wavelength[i] - filter.wavelength[i-1])
            tr = filter.transmission[i-1] * lw + filter.transmission[i] * (1 - lw)
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
