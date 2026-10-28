# DynamicProps - a property bag that supports dot-notation access (like MATLAB structs / Python SimpleNamespace)
#
# Usage:
#   p = DynamicProps()
#   p.alpha = 2.0          # set
#   p.alpha                # get => 2.0
#   haskey(p, :alpha)      # => true
#   isempty(p)             # => true/false
mutable struct DynamicProps
    _data::Dict{Symbol, Any}
end
DynamicProps() = DynamicProps(Dict{Symbol, Any}())

function Base.getproperty(dp::DynamicProps, s::Symbol)
    s === :_data && return getfield(dp, :_data)
    d = getfield(dp, :_data)
    haskey(d, s) || error("DynamicProps has no property '$s'")
    return d[s]
end

function Base.setproperty!(dp::DynamicProps, s::Symbol, v)
    s === :_data && return setfield!(dp, :_data, v)
    getfield(dp, :_data)[s] = v
end

Base.haskey(dp::DynamicProps, s::Symbol) = haskey(getfield(dp, :_data), s)
Base.isempty(dp::DynamicProps) = isempty(getfield(dp, :_data))
