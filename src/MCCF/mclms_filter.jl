function mclms_filter(
  X_r::AbstractVector,
  X_r_prime::AbstractVector,
  mclms_len::Integer,
)
  length(X_r) == length(X_r_prime) ||
    throw(DimensionMismatch("input vectors must have equal lengths"))
  mclms_len > 0 || throw(ArgumentError("mclms_len must be positive"))
  n = length(X_r)
  denom = sum(abs2, X_r)
  μ = denom == 0 ? 0.0 : 1 / denom
  T = promote_type(eltype(X_r), eltype(X_r_prime), ComplexF64)
  weights = zeros(T, mclms_len)
  xvec = zeros(T, mclms_len)
  for j = 1:n
    fill!(xvec, 0)
    lo = max(1, j-mclms_len+1)
    count = j-lo+1
    xvec[(mclms_len - count + 1):mclms_len] .= X_r[lo:j]
    y = sum(weights .* xvec)
    err = X_r_prime[j] - y
    weights .= 0.98 .* weights .+ μ .* err .* conj.(xvec)
  end
  # MATLAB uses filter(fliplr(W), 1, X_r) after training.
  b = reverse(weights)
  out = zeros(T, n)
  for j = 1:n
    upto = min(mclms_len, j)
    for k = 1:upto
      out[j] += b[k] * X_r[j - k + 1]
    end
  end
  return out
end
