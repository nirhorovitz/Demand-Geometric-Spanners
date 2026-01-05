module IOUtils

using ..CoreTypes
using JLD2
using FileIO

export save_results, load_results

"""
    save_results(filename::String, instance::SpannerInstance, results::Vector{SpannerResult})

Saves the instance and results to a JLD2 file.
"""
function save_results(filename::String, instance::SpannerInstance, results::Vector{SpannerResult})
    # JLD2 handles Julia structs efficiently
    save(filename, Dict("instance" => instance, "results" => results))
end

"""
    load_results(filename::String)

Loads results from a file.
"""
function load_results(filename::String)
    data = load(filename)
    return data["instance"], data["results"]
end

end

