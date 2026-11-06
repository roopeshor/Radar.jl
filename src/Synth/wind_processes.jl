function convective_plumes(
  h::AbstractVector;
  plume_count::Integer = 4,
  plume_max_alt::Real = 6.0,
)
  z = Float64.(h)
  out = zeros(length(z))
  for _ = 1:plume_count
    center = 0.8 + rand()*(plume_max_alt-0.8)
    peak = (rand()-0.3)*2.5
    width = 0.2 + rand()*0.4
    out .+= peak .* exp.(-((z .- center) ./ width) .^ 2)
  end
  return out
end

function gravity_waves(
  h::AbstractVector;
  gw_nmodes::Integer = 20,
  gw_amp::Real = 2.5,
  gw_max_lambda::Real = 4.0,
  gw_damp_scale::Real = 12.0,
  gw_cf::Real = 0.05,
  lw_Hs::Real = 7.0,
)
  z = Float64.(h)
  n = length(z)
  m = collect(range(2π/gw_max_lambda, 2π/0.4, length = gw_nmodes))
  phase_u = 2π .* rand(gw_nmodes)
  phase_v = 2π .* rand(gw_nmodes)
  energy = m .^ (-1.5)
  energy ./= sum(energy)
  uw = zeros(n)
  vw = zeros(n)
  for i = 1:gw_nmodes
    uw .+= sqrt(energy[i]) .* cos.(m[i] .* z .+ phase_u[i])
    vw .+= sqrt(energy[i]) .* sin.(m[i] .* z .+ phase_v[i])
  end
  amp = exp.(z ./ (2lw_Hs)) .* exp.(-z ./ gw_damp_scale)
  su = std(uw)
  sv = std(vw)
  su > 0 && (uw ./= su)
  sv > 0 && (vw ./= sv)
  u = gw_amp .* uw .* amp
  v = gw_amp .* vw .* amp
  return (u, v, -gw_cf .* u)
end

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

function orographic_updraft(
  h::AbstractVector,
  u_mean::AbstractVector;
  ou_alpha_slope::Real = atan(3/20),
  ou_H_lift::Real = 2.0,
)
  return max.(0.0, abs.(u_mean)) .* tan(ou_alpha_slope) .* exp.(-h ./ ou_H_lift)
end

function somali_jet(
  h::AbstractVector;
  sj_w::Real = 12.0,
  sj_uv::Real = 14.0,
  sj_center::Real = 1.5,
  sj_s::Real = 0.8,
)
  shape = exp.(-((h .- sj_center) ./ sj_s) .^ 2)
  return (sj_w .* shape, sj_uv .* shape)
end

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

function synoptic_bg(
  h::AbstractVector;
  bg_z0_m::Real = 0.05,
  bg_ustar::Real = 0.3,
  bg_kappa::Real = 0.4,
  bg_blh::Real = 1.5,
)
  z_m = h .* 1000
  return (bg_ustar/bg_kappa) .* log.((z_m .+ bg_z0_m) ./ bg_z0_m) .*
         exp.(-h ./ bg_blh)
end

function tropical_easterly_jet(
  h::AbstractVector;
  tej_u::Real = 10.0,
  tej_v::Real = 3.0,
  tej_alt::Real = 15.0,
  tej_vthick::Real = 2.5,
)
  shape = exp.(-((h .- tej_alt) ./ tej_vthick) .^ 2)
  return (tej_u .* shape, tej_v .* shape)
end
