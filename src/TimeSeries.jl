using Reexport
@reexport module TimeSeries

abstract type GenericTimeSeries end

struct SingleTimeSeries{T} <: GenericTimeSeries
    ts::Vector{T}
end

export SingleTimeSeries
end
