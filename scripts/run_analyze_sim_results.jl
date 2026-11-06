using Statistics
using Plots

function read_tsv(path)
  lines   = readlines(path)
  headers = split(first(lines), '\t')
  rows    = [split(line, '\t') for line in lines[2:end] if !isempty(strip(line))]
  return headers, rows
end

function main(folder = "outputs/tester")
  mh, mrows = read_tsv(joinpath(folder, "mse.tsv"))
  ch, crows = read_tsv(joinpath(folder, "correlation.tsv"))
  mi = Dict(k=>findfirst(==(k), mh) for k in mh)
  ci = Dict(k=>findfirst(==(k), ch) for k in ch)

  methods     = sort(unique(row[mi["method"]] for row in mrows))
  beams       = sort(unique(row[mi["direction"]] for row in mrows))
  mean_mse    = fill(NaN, length(methods), length(beams))
  median_mse  = similar(mean_mse)
  mean_corr   = similar(mean_mse)
  median_corr = similar(mean_mse)

  for (i, method) in enumerate(methods), (j, beam) in enumerate(beams)
    vals = [
      parse(Float64, row[mi["MSE"]]) for
      row in mrows if row[mi["method"]]==method && row[mi["direction"]]==beam
    ]
    cors = [
      parse(Float64, row[ci["Correlation"]]) for
      row in crows if row[ci["method"]]==method && row[ci["direction"]]==beam
    ]
    !isempty(vals) &&
      (mean_mse[i, j] = mean(vals); median_mse[i, j] = median(vals))
    !isempty(cors) &&
      (mean_corr[i, j] = mean(cors); median_corr[i, j] = median(cors))
  end
  display((
    methods     = methods,
    beams       = beams,
    mean_mse    = mean_mse,
    median_mse  = median_mse,
    mean_corr   = mean_corr,
    median_corr = median_corr,
  ))
  p1 = bar(
    methods,
    mean(mean_corr; dims = 2)[:];
    ylabel = "Mean correlation",
    title = "Mean correlation by method",
    legend = false,
  )
  p2 = bar(
    methods,
    mean(mean_mse; dims = 2)[:];
    yscale = :log10,
    ylabel = "Mean MSE",
    title = "Mean MSE by method",
    legend = false,
  )
  display(plot(p1, p2; layout = (1, 2), size = (1000, 450)))
  return (; methods, beams, mean_mse, median_mse, mean_corr, median_corr)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
  main()
end
