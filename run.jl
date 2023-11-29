include("src/main.jl")
using Plots
using ProgressMeter
using DelimitedFiles
import UnitfulAstro: ly

FILTER_TYPE = :photon

ser = read_series_filterdata(
    FilterFolder("data/Filters/", FILTER_TYPE),
    "data/data_13dqy_formatted_for_package.txt",
    unit=Flux(160e6ly))

dates = readdlm("data/13dqy_int_dates.txt") |> vec

fser = fit(ser, BlackBodyModel(), dates)
plot(fser)
