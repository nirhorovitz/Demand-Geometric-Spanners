module SpannerComparison

using Reexport

include("Core.jl")
@reexport using .CoreTypes

include("Generators.jl")
@reexport using .Generators

include("Algorithms.jl")
@reexport using .Algorithms

include("Analysis.jl")
@reexport using .Analysis

include("Visualization.jl")
@reexport using .Visualization

include("IO.jl")
@reexport using .IOUtils

end

