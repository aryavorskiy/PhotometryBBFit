using Downloads, Logging

# Getting filter info

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

get_filter(a::Any, ::String) = error("`$a` is not a valid filter provider")
get_filter(f::Function, tag::String) = f(tag)

struct SVO2Website end
function get_filter(::SVO2Website, id::String)
    mode=:photon
    file = tempname()
    Downloads.download("http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=$id", file)
    filter = read_filter(file, mode, id)
    rm(file)
    return filter
end

struct FilterFolder
    path::String
    type::Symbol
end
function get_filter(fl::FilterFolder, id::String)
    candidates = []
    for name in readdir(fl.path)
        mt = match(r"^([\w\.]*)?(\.(rtf|txt|dat|csv))$", name)
        mt === nothing && continue
        if mt[1] == id
            push!(candidates, name)
        end
    end
    if isempty(candidates)
        error("Filter `$id` not found in folder $(fl.path)")
    elseif length(candidates) ≥ 2
        error("""Filter `$id` is ambiguous in folder `$(fl.path)`:
        $(join(candidates, ", "))""")
    end
    return read_filter(joinpath(fl.path, only(candidates)), fl.type, id)
end

# Reading data series

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

function read_series_filterdata(callback, file; unit=Luminosity())
    tags, n_sers = read_series(file)
    filters = [get_filter(callback, tag) for tag in tags]
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
