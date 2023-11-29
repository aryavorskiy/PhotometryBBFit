using PhotometryFit
using Plots
using DelimitedFiles
import UnitfulAstro: ly

ser = read_series_filterdata(
    FilterFolder("data/Filters/", :photon),
    "data/data_13dqy_formatted_for_package.txt",
    unit=Flux(dist=160e6ly))

dates = readdlm("data/13dqy_int_dates.txt") |> vec

fser = fit(ser, BlackBodyModel(), dates)
scatter(fser)
