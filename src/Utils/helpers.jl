"MATLAB-compatible FFT shift along one or more dimensions."
function fftshift(x::AbstractArray, dims = 1:ndims(x))
  ds = dims isa Integer ? (dims,) : Tuple(dims)
  shifts = ntuple(d -> d in ds ? fld(size(x, d), 2) : 0, ndims(x))
  return circshift(x, shifts)
end

"Inverse of `fftshift`, including odd-length dimensions."
function ifftshift(x::AbstractArray, dims = 1:ndims(x))
  ds = dims isa Integer ? (dims,) : Tuple(dims)
  shifts = ntuple(d -> d in ds ? cld(size(x, d), 2) : 0, ndims(x))
  return circshift(x, shifts)
end

function hanning_window(n::Integer)
  n <= 0 && return Float64[]
  n == 1 && return ones(1)
  return 0.5 .- 0.5 .* cos.(2pi .* (0:(n - 1)) ./ (n-1))
end

"Centered moving average with truncated windows at array boundaries."
function moving_mean(x::AbstractVector, width::Integer)
  width <= 1 && return copy(x)
  n = length(x)
  y = similar(float.(x))
  left = fld(width - 1, 2)
  right = width - left - 1
  for i = 1:n
    lo = max(1, i-left)
    hi = min(n, i+right)
    y[i] = mean(@view x[lo:hi])
  end
  return y
end
function moving_mean(x::AbstractMatrix, width::Integer; dims::Integer = 2)
  dims == 2 || throw(ArgumentError("moving_mean currently supports dims=2"))
  y = similar(float.(x))
  for r in axes(x, 1)
    y[r, :] = moving_mean(@view(x[r, :]), width)
  end
  return y
end

"Gaussian smoothing with a normalized finite kernel and replicated edge values."
function gaussian_smooth(x::AbstractVector, width::Integer)
  width <= 1 && return copy(x)
  w = max(3, width)
  w = iseven(w) ? w + 1 : w
  σ = max(0.5, w / 6)
  radius = fld(w, 2)
  kernel = exp.(-(((-radius):radius) .^ 2) ./ (2σ^2))
  kernel ./= sum(kernel)
  y = similar(float.(x))
  n = length(x)
  for i = 1:n
    acc = zero(eltype(y))
    for k = (-radius):radius
      j = clamp(i+k, 1, n)
      acc += kernel[k + radius + 1] * x[j]
    end
    y[i] = acc
  end
  return y
end
function gaussian_smooth(x::AbstractMatrix, width::Integer; dims::Integer = 2)
  dims == 2 || throw(ArgumentError("gaussian_smooth currently supports dims=2"))
  y = similar(float.(x))
  for r in axes(x, 1)
    y[r, :] = gaussian_smooth(@view(x[r, :]), width)
  end
  return y
end

"Indices of local maxima, excluding endpoints as MATLAB `findpeaks` does."
function local_peaks(x::AbstractVector)
  n = length(x)
  n < 3 && return Int[]
  idx = Int[]
  i = 2
  while i < n
    if x[i] > x[i - 1]
      # Collapse a flat-topped peak to its middle sample.
      j = i
      while j < n && x[j] == x[j + 1]
        j += 1
      end
      if j < n && x[j] > x[j + 1]
        push!(idx, fld(i+j, 2))
      end
      i = j + 1
    else
      i += 1
    end
  end
  return idx
end

"Linearly map values from one interval to another."
function rescale_linear(
  x,
  out_min::Real,
  out_max::Real,
  in_min::Real,
  in_max::Real,
)
  in_max == in_min && return fill((out_min + out_max) / 2, size(x))
  return out_min .+ (x .- in_min) .* ((out_max - out_min) / (in_max - in_min))
end

function safe_corrcoef(x::AbstractVector, y::AbstractVector)
  length(x) == length(y) ||
    throw(DimensionMismatch("vectors must have equal length"))
  mask = isfinite.(x) .& isfinite.(y)
  count(mask) < 2 && return NaN
  a = x[mask]
  b = y[mask]
  (std(a) == 0 || std(b) == 0) && return NaN
  return cor(a, b)
end

"""
    fft_n(x, n)

Calculate an padded FFT of length `n`, following MATLAB's fft(x, n) behavior.

- If n > length(x), zero-pad the input.
- If n < length(x), truncate the input.
- If n == length(x), perform the ordinary FFT.
"""
function fft_n(x::AbstractVector, n::Integer)
  n > 0 || throw(ArgumentError("FFT length must be positive"))

  len = length(x)

  if n == len
    return fft(x)
  elseif n < len
    return fft(view(x, 1:n))
  else
    padded = zeros(eltype(x), n)
    padded[1:len] .= x
    return fft(padded)
  end
end

function fft_n(x::AbstractArray, n::Integer, dim::Integer)
  n > 0 || throw(ArgumentError("FFT length must be positive"))
  1 <= dim <= ndims(x) || throw(ArgumentError("Invalid FFT dimension"))

  old_len = size(x, dim)

  if n == old_len
    return fft(x, dim)
  end

  new_dims = ntuple(d -> d == dim ? n : size(x, d), ndims(x))

  padded = zeros(eltype(x), new_dims)

  copy_len = min(n, old_len)
  inds = ntuple(d -> d == dim ? (1:copy_len) : Colon(), ndims(x))

  padded[inds...] .= view(x, inds...)

  return fft(padded, dim)
end


function gradient(x::AbstractVector, h::Real)
  n = length(x)
  n < 2 && return zeros(Float64, n)
  g = similar(Float64.(x))
  g[1] = (x[2]-x[1])/h
  g[end] = (x[end]-x[end - 1])/h
  for i = 2:(n - 1)
    g[i] = (x[i + 1]-x[i - 1])/(2h)
  end
  return g
end
