function stratospheric_qbo(
  h::AbstractVector;
  qbo_lambda::Real = 20.0,
  qbo_u::Real = 8.0,
  qbo_v::Real = 3.0,
  qbo_alt::Real = 17.0,
)
  envelope = 0.5 .* (1 .+ tanh.((h .- qbo_alt) ./ 2))
  phase = 2π .* (h .- qbo_alt) ./ qbo_lambda
  return (qbo_u .* sin.(phase) .* envelope, qbo_v .* cos.(phase) .* envelope)
end
