module Visualization

using ..CoreTypes
using Plots
using Graphs
using SimpleWeightedGraphs

export visualize_results

"""
    visualize_results(instance::SpannerInstance, results::Vector{SpannerResult})

Plots the resulting graphs and statistics.
"""
function visualize_results(instance::SpannerInstance, results::Vector{SpannerResult})
    # Extract points for plotting
    xs = [p[1] for p in instance.points]
    ys = [p[2] for p in instance.points]
    
    plots = []
    
    for res in results
        g = res.graph
        
        # Build comprehensive stats text
        s = res.stats
        title_text = "$(res.algorithm_name)\n" *
                     "Time: $(round(res.runtime_seconds, digits=3))s\n" *
                     "Valid: $(get(s, :is_valid_spanner, "?")) (t=$(get(s, :target_t, "?")))\n" *
                     "Weight: $(round(s[:total_weight], digits=1)) (x$(round(get(s, :weight_ratio, 0), digits=2)))\n" *
                     "Edges: $(s[:num_edges]) | Pts: $(s[:num_nodes])\n" *
                     "Deg: Max $(s[:max_degree]) / Avg $(round(s[:avg_degree], digits=1))"
        
        p = plot(xs, ys, seriestype=:scatter, markersize=2, legend=false, title=title_text, titlefontsize=10, aspect_ratio=:equal)
        
        # Draw edges
        for e in edges(g)
            u, v = src(e), dst(e)
            plot!(p, [xs[u], xs[v]], [ys[u], ys[v]], color=:gray, alpha=0.5)
        end
        
        push!(plots, p)
    end
    
    # Combine plots
    final_plot = plot(plots..., layout=(1, length(results)), size=(400*length(results), 600), margin=5Plots.mm)
    return final_plot
end

end

