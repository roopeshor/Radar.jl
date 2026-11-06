function kolmogorov_turbulence(
  h::AbstractVector,
  u_mean::AbstractVector,
  v_mean::AbstractVector;
  kt_vaf::Real = 0.25,
  kt_bf::Real = 1.0,
  kt_maxT::Real = 2.5,
  kt_Ric::Real = 0.25,
  kt_tlpeak::Real = 3.0,
  kt_vthick::Real = 2.5,
  turb_interm_sigma::Real = 0.5,
  lw_N_bv::Real = 0.012,
)
  z = Float64.(h)
  n = length(z)
  n >= 2 || throw(ArgumentError("At least two altitude bins are required"))
  dz_m = (z[2]-z[1])*1000
  du = gradient(Float64.(u_mean), dz_m)
  dv = gradient(Float64.(v_mean), dz_m)
  shear2 = du .^ 2 .+ dv .^ 2 .+ 1e-8
  ri = lw_N_bv^2 ./ shear2
  boost =
    kt_bf .+
    kt_maxT .* (ri .< kt_Ric) .* exp.(-((z .- kt_tlpeak) ./ kt_vthick) .^ 2)
  k =
    vcat(collect(0:floor(Int, n / 2)), collect((-floor(Int, (n - 1) / 2)):-1)) ./
    (n*(z[2]-z[1]))
  kmin = 1/(n*(z[2]-z[1]))
  kabs = max.(abs.(k), kmin)
  transfer = kabs .^ (-5/6)
  transfer[1] = 0
  rawu = real.(ifft(fft(randn(n)) .* transfer))
  rawv = real.(ifft(fft(randn(n)) .* transfer))
  raww = real.(ifft(fft(randn(n)) .* transfer))
  normu = (rawu .- mean(rawu)) ./ max(std(rawu), eps())
  normv = (rawv .- mean(rawv)) ./ max(std(rawv), eps())
  normw = (raww .- mean(raww)) ./ max(std(raww), eps())
  μ = -0.5*turb_interm_sigma^2
  interm = exp.(μ .+ turb_interm_sigma .* normu)
  base = 0.15 .+ 0.4 .* sqrt.(shear2) ./ maximum(sqrt.(shear2))
  sigma = base .* boost .* interm
  return (normu .* sigma, normv .* sigma, normw .* (sigma .* kt_vaf))
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
