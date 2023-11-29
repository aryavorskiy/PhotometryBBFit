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
        wl, ts = parse.(Float64, split(line, r"\s+|\s*,\s*"))
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
    mode = :photon  # TODO: download real data from the site
    file = tempname()
    Downloads.download("http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=$id", file)
    filter = read_filter(file, mode, id)
    rm(file)
    return filter
end

struct FilterFolder
    path::String
    default::Symbol
    dict::Dict{String, Tuple{String, Symbol}}
end
FilterFolder(path::String, default::Symbol=:photon) =
    FilterFolder(path, default, Dict{String, Tuple{String, Symbol}}())
function get_filter(fl::FilterFolder, id::String)
    candidates = []
    for name in readdir(fl.path)
        new_id, mode = get(fl.dict, id, (id, fl.default))
        mode = mode === :default ? fl.default : mode
        mt = match(r"^([\w\.]*)?(\.(rtf|txt|dat|csv))$", name)
        mt === nothing && continue
        if mt[1] == new_id
            push!(candidates, (name, mode, new_id))
        end
    end
    if isempty(candidates)
        return nothing
    elseif length(candidates) ≥ 2
        error("""Filter `$id` is ambiguous in folder `$(fl.path)`:
        $(join(candidates, "\n"))""")
    end
    return read_filter(joinpath(fl.path, only(candidates)[1]), only(candidates)[2], only(candidates)[3])
end

function ConversionRules(pairs::Pair{String}...)
    d = Dict{String, Tuple{String, Symbol}}()
    for (l, r) in pairs
        if r isa Tuple{String, Symbol}
            d[l] = r
        elseif r isa Tuple{Symbol, String}
            d[l] = (r[2], r[1])
        elseif r isa String
            d[l] = (r, :default)
        elseif r isa Symbol
            d[l] = (l, r)
        else
            error("Invalid conversion rule $r")
        end
    end
    return d
end

# Reading data series

function read_photometry_data(provider, file; unit=Luminosity())
    filters_and_series = Dict{Filter{Float64}, Series{Float64}}()
    unresolved_tags = String[]
    for line in eachline(file)
        startswith(line, '#') && continue
        time_s, val_s, err_s, tag = split(line, r"\s+|\s*,\s*")
        time, val, err = parse.(Float64, (time_s, val_s, err_s))

        # resolve filter
        filter = get_filter(provider, String(tag))
        if filter === nothing
            tag in unresolved_tags || push!(unresolved_tags, tag)
            continue
        end
        if !(filter in keys(filters_and_series))
            filters_and_series[filter] = Series()
        end

        insert_measurement!(filters_and_series[filter], time, val, err)
    end
    filters = collect(keys(filters_and_series))
    sers = [to_luminosity(unit, ser) for ser in values(filters_and_series)]
    sp = sortperm(filters, by=lambda_eff)
    return PhotometryData(filters[sp], sers[sp])
end
