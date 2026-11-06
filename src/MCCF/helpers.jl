function _conv1d(x::AbstractVector, y::AbstractVector)
  out = zeros(promote_type(eltype(x), eltype(y)), length(x)+length(y)-1)
  for i in eachindex(x), j in eachindex(y)
    out[i + j - 1] += x[i] * y[j]
  end
  return out
end

function _normalise_spectra(spectra::AbstractMatrix, axis)
  out = Float64.(real.(spectra))
  if axis === :all || axis == 0 || axis == "all"
    mx = maximum(out)
    return mx == 0 ? out : out ./ mx
  end
  mx = maximum(out; dims = axis)
  mx[mx .== 0] .= 1
  return out ./ mx
end

function _remove_dc_bins!(f_sig, noise_val)
  nfft = size(f_sig, ndims(f_sig))
  mask = [1, 2, 3, nfft, nfft-1]
  if ndims(f_sig) == 1
    f_sig[mask] .= noise_val .* cis.(2π .* rand(length(mask)))
  else
    for r in axes(f_sig, 1)
      f_sig[r, mask] .= noise_val[r] .* cis.(2π .* rand(length(mask)))
    end
  end
  return f_sig
end
