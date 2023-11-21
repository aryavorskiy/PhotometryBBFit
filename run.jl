include("src/main.jl")
using Plots
import UnitfulAstro: ly

ser = read_series_filterdata("data_13dqy_formatted_for_package.txt", unit=Flux(16e10ly)) do filter
    if isfile("Filters/$filter.txt")
        return read_filter("Filters/$filter.txt")
    elseif isfile("Filters/$filter.rtf")
        return read_filter("Filters/$filter.rtf")
    else
        return nothing
    end
end

sp = ser(5)
e = 10.0^-18
Rs = e:e:100e
Ts = 4000:1000:8000
data = [chi2(sp, (R, T)) for R in Rs, T in Ts]
heatmap(Rs, Ts, -data)
