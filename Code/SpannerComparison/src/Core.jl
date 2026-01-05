module CoreTypes

using StaticArrays
using Graphs
using SimpleWeightedGraphs

export Point2D, SpannerInstance, SpannerResult, AbstractSpannerAlgorithm, MatrixWFunc

const Point2D = SVector{2, Float64}

"""
    MatrixWFunc

A callable struct that wraps a weight matrix to be used as a function.
This ensures serializability with JLD2 compared to anonymous functions.
"""
struct MatrixWFunc <: Function
    W::Matrix{Float64}
end

(m::MatrixWFunc)(i::Int, j::Int) = m.W[i, j]

"""
    SpannerInstance

Holds the problem instance data.
- `points`: Vector of 2D points.
- `w_func`: A function that takes two indices (i, j) and returns a weight w(p_i, p_j).
- `t`: The stretch factor parameter.
"""
struct SpannerInstance{F<:Function}
    points::Vector{Point2D}
    w_func::F
    t::Float64
end

"""
    AbstractSpannerAlgorithm

Abstract type for all spanner construction algorithms.
"""
abstract type AbstractSpannerAlgorithm end

"""
    SpannerResult

Holds the result of a spanner algorithm.
- `algorithm_name`: Name of the algorithm used.
- `graph`: The resulting graph (weighted).
- `runtime_seconds`: Time taken to compute the graph.
- `stats`: A dictionary of computed statistics.
"""
struct SpannerResult
    algorithm_name::String
    graph::SimpleWeightedGraph{Int, Float64}
    runtime_seconds::Float64
    stats::Dict{Symbol, Any}
end

end

