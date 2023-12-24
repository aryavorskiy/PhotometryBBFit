using Downloads, Logging

# Getting filter info

"""
    read_filter(file)

Reads filter info from `file` in tab-separated format. Returns a `Filter`.
"""
function read_filter(file, meta=FilterMetadata())
    wavelength = Float64[]
    transparency = Float64[]
    for line in eachline(file)
        wl, ts = parse.(Float64, split(line, r"\s+|\s*,\s*"))
        push!(wavelength, wl)
        push!(transparency, ts)
    end
    length(wavelength) == 0 && @warn "Empty filter data in filename `$file`"
    return Filter{Float64}(wavelength, transparency, meta)
end

get_filter(a::Any, ::String) = error("`$a` is not a valid filter provider")
get_filter(f::Function, tag::String) = f(tag)

# struct SVO2Website end
# function get_filter(::SVO2Website, id::String)
#     mode = :photon  # TODO: download real data from the site
#     file = tempname()
#     Downloads.download("http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=$id", file)
#     filter = read_filter(file, mode, id)
#     rm(file)
#     return filter
# end

struct FilterFolder
    path::String
    rules::Dict{String, Tuple{String, FilterMetadata}}
    default_meta::FilterMetadata
end
function FilterInfo(path=nothing, band::Symbol=:none; kw...)
    return path, FilterMetadata(;band=band, kw...)
end
function FilterFolder(path::String, rules::Pair{String, <:Tuple}...; kw...)
    rs = Dict{String, Tuple{String, FilterMetadata}}()
    for (tag, (path, meta)) in rules
        if path === nothing
            path = tag
        end
        new_meta = FilterMetadata(;
            kw..., id=tag, family=meta.family, system=meta.system, band=meta.band, mode=meta.mode)
        rs[tag] = (path, new_meta)
    end
    FilterFolder(path, rs, FilterMetadata(;kw...))
end
function get_filter(fl::FilterFolder, tag::String)
    candidates = []
    filename, meta = get(fl.rules, tag, (tag, fl.default_meta))
    for name in readdir(fl.path)
        mt = match(r"^([\w\.]*)?(\.(rtf|txt|dat|csv))$", name)
        mt === nothing && continue
        if mt[1] == filename
            push!(candidates, name)
        end
    end
    if isempty(candidates)
        return nothing
    elseif length(candidates) ≥ 2
        error("""Filter '$tag' is ambiguous in folder `$(fl.path)`:
        $(join(candidates, "\n"))""")
    end
    return read_filter(joinpath(fl.path, only(candidates)), meta)
end

# Reading data series

function read_photometry_data(provider, file; unit=Luminosity())
    filters_and_series = Dict{Filter{Float64}, Series{Float64}}()
    unresolved_filters = String[]
    for line in eachline(file)
        startswith(line, '#') && continue
        time_s, val_s, err_s, tag = split(line, r"\s+|\s*,\s*")
        time, val, err = parse.(Float64, (time_s, val_s, err_s))

        # resolve filter
        filter = get_filter(provider, String(tag))
        if filter === nothing
            tag in unresolved_filters || push!(unresolved_filters, tag)
            continue
        end
        if !(filter in keys(filters_and_series))
            filters_and_series[filter] = Series()
        end

        insert_measurement!(filters_and_series[filter], time, val, err)
    end
    isempty(unresolved_filters) || @warn unresolved_filters
    filters = collect(keys(filters_and_series))
    sers = [to_luminosity(unit, filter, ser) for (filter, ser) in filters_and_series]
    return PhotometryData{Float64}(filters, sers)
end
