# The instrument format is little-endian, independent of host byte order.
@inline function _readle(::Type{Int16}, b::AbstractVector{UInt8}, i::Int)
  u = UInt16(b[i]) | (UInt16(b[i + 1]) << 8)
  return reinterpret(Int16, u)
end
@inline function _readle(::Type{UInt32}, b::AbstractVector{UInt8}, i::Int)
  u =
    UInt32(b[i]) | (UInt32(b[i + 1]) << 8) | (UInt32(b[i + 2]) << 16) |
    (UInt32(b[i + 3]) << 24)
  return u
end
@inline _readle(::Type{Int32}, b::AbstractVector{UInt8}, i::Int) =
  reinterpret(Int32, _readle(UInt32, b, i))
@inline function _readle(::Type{Float32}, b::AbstractVector{UInt8}, i::Int)
  return reinterpret(Float32, _readle(UInt32, b, i))
end
@inline _readle(::Type{Int8}, b::AbstractVector{UInt8}, i::Int) =
  reinterpret(Int8, b[i])

@inline function _putle!(b::Vector{UInt8}, i::Int, x::Int16)
  u = reinterpret(UInt16, x)
  b[i] = UInt8(u & 0xff)
  b[i + 1] = UInt8((u >> 8) & 0xff)
end
@inline function _putle!(b::Vector{UInt8}, i::Int, x::UInt32)
  b[i] = UInt8(x & 0xff)
  b[i + 1] = UInt8((x >> 8) & 0xff)
  b[i + 2] = UInt8((x >> 16) & 0xff)
  b[i + 3] = UInt8((x >> 24) & 0xff)
end
@inline function _putle!(b::Vector{UInt8}, i::Int, x::Int32)
  _putle!(b, i, reinterpret(UInt32, x))
end
@inline function _putle!(b::Vector{UInt8}, i::Int, x::Float32)
  _putle!(b, i, reinterpret(UInt32, x))
end
@inline function _putle!(b::Vector{UInt8}, i::Int, x::Int8)
  b[i] = reinterpret(UInt8, x)
end

function _put_vector!(b, start, values, ::Type{T}) where {T}
  width = sizeof(T)
  for (j, value) in enumerate(values)
    _putle!(b, start + (j-1)*width, convert(T, value))
  end
  return b
end

function _fixed_bytes(value, n::Integer)
  bytes = value isa AbstractString ? collect(codeunits(value)) : UInt8.(value)
  out = zeros(UInt8, n)
  copyto!(out, 1, bytes, 1, min(length(bytes), n))
  return out
end

@inline _round_away_int32(x::Real) =
  x >= 0 ? floor(Int32, x + 0.5) : ceil(Int32, x - 0.5)
