"Hildebrand–Sekhon noise estimator, operating row-wise on power spectra."
function HS_noise_estimate(spectra::AbstractMatrix, in_coh::Real)
  nr, nfft = size(spectra)
  noise = zeros(Float64, nr)
  m = max(Float64(in_coh), 1.0)
  for r = 1:nr
    a = sort(Float64.(real.(spectra[r, :])))
    pn = cumsum(a) ./ (1:nfft)
    second = cumsum(abs2.(a)) ./ (1:nfft)
    qn = second .- pn .^ 2
    valid = findall((qn .> 0) .& ((pn .^ 2 ./ (m .* qn)) .> 1))
    k = isempty(valid) ? max(1, round(Int, 0.1*nfft)) : last(valid)
    noise[r] = pn[k]
  end
  return noise
end
