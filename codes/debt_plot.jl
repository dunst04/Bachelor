using CSV   
using DataFrames
using Plots
using Dates

# 1. Load the CSV file into a DataFrame
csv_path = joinpath(@__DIR__, "debt_to_gdp.csv")
df = CSV.read(csv_path, DataFrame)
dropmissing!(df)

# 2. Parse the dates so the x-axis acts as a proper time series
# Based on the snippet, dates look like "01/01/1998" (mm/dd/yyyy)
dates = Date.(df.observation_date, "mm/dd/yyyy")

# 3. Create the plot with the specific colors and labels
p = plot(yformatter = y -> "$(round(Int, y))%")

# US data: green, label "USA"
plot!(p, dates, df.us, label="USA", color=:green, linewidth=2)

# Euro area data: blue, label "euro area"
# Note: Since the column name has a space, we use df[!, "euro area"]
plot!(p, dates, df[!, "euro area"], label="euro area", color=:blue, linewidth=2)

# Japan data: red, label "Japan"
plot!(p, dates, df.japan, label="Japan", color=:red, linewidth=2)

# Display the final plot
display(p)