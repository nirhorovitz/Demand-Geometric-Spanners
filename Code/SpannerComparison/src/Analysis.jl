module Analysis

using ..CoreTypes
using Graphs
using SimpleWeightedGraphs
using LinearAlgebra
using Statistics
using Distances

export compute_stats

"""
    compute_stats(instance::SpannerInstance, result::SpannerResult)

Computes statistics for the generated graph and updates the result.stats dictionary.
Uses parallel processing for spanner verification.
"""
function compute_stats(instance::SpannerInstance, result::SpannerResult)
    g = result.graph
    points = instance.points
    n = length(points)
    
    # 1. Basic Graph Stats
    num_edges = ne(g)
    num_nodes = nv(g)
    degrees = degree(g)
    max_deg = isempty(degrees) ? 0 : maximum(degrees)
    avg_deg = isempty(degrees) ? 0.0 : mean(degrees)
    
    # 2. Total Weight
    total_weight = sum(weight(e) for e in edges(g))
    
    # 3. MST Weight
    mst_weight = compute_mst_weight(points)
    weight_multiplier = mst_weight > 0 ? total_weight / mst_weight : Inf
    
    # 4. Spanner Verification (Parallel)
    max_stretch_found, is_valid_spanner = verify_spanner_properties(instance, g)
    
    # Update stats
    result.stats[:total_weight] = total_weight
    result.stats[:mst_weight] = mst_weight
    result.stats[:weight_ratio] = weight_multiplier
    result.stats[:num_edges] = num_edges
    result.stats[:num_nodes] = num_nodes
    result.stats[:max_degree] = max_deg
    result.stats[:avg_degree] = avg_deg
    result.stats[:target_t] = instance.t
    result.stats[:max_stretch_found] = max_stretch_found
    result.stats[:is_valid_spanner] = is_valid_spanner
    
    return result
end

function verify_spanner_properties(instance::SpannerInstance, g::SimpleWeightedGraph)
    n = length(instance.points)
    points = instance.points
    w_func = instance.w_func
    t = instance.t
    
    # Shared variables with lock
    max_stretch_found = Ref(0.0)
    is_valid_spanner = Ref(true)
    lock_obj = ReentrantLock()
    
    # Parallelize outer loop (source nodes)
    Threads.@threads for i in 1:n
        # Calculate shortest paths from node i to all other nodes
        dists = dijkstra_shortest_paths(g, i).dists
        
        local_max = 0.0
        local_valid = true
        
        for j in (i+1):n
            d_graph = dists[j]
            d_euclid = norm(points[i] - points[j])
            
            if d_euclid > 1e-10
                # Condition: w(p,q) * d_G(p,q) <= t * ||p-q||
                # Effective Stretch: (w(p,q) * d_G(p,q)) / ||p-q||
                w_val = w_func(i, j)
                
                # If d_graph is Inf (not connected), stretch is Inf
                if d_graph == Inf
                    stretch_val = Inf
                else
                    stretch_val = (w_val * d_graph) / d_euclid
                end
                
                if stretch_val > local_max
                    local_max = stretch_val
                end
                
                # Check validity with small tolerance for float errors
                if stretch_val > t + 1e-7
                    local_valid = false
                end
            end
        end
        
        # Update shared state
        lock(lock_obj) do
            if local_max > max_stretch_found[]
                max_stretch_found[] = local_max
            end
            if !local_valid
                is_valid_spanner[] = false
            end
        end
    end
    
    return max_stretch_found[], is_valid_spanner[]
end

function compute_mst_weight(points)
    n = length(points)
    if n == 0 return 0.0 end
    
    min_dists = fill(Inf, n)
    in_mst = fill(false, n)
    min_dists[1] = 0.0
    total_w = 0.0
    
    for _ in 1:n
        u = -1
        min_val = Inf
        
        # Find min distance node not in MST
        for i in 1:n
            if !in_mst[i] && min_dists[i] < min_val
                min_val = min_dists[i]
                u = i
            end
        end
        
        if u == -1 break end
        
        in_mst[u] = true
        total_w += min_val
        
        # Update neighbors (implicit complete graph)
        # Using squared distance might be faster for comparison, but we need actual weight sum.
        # Threads.@threads not useful for inner loop usually unless N is very huge.
        for v in 1:n
            if !in_mst[v]
                d = norm(points[u] - points[v])
                if d < min_dists[v]
                    min_dists[v] = d
                end
            end
        end
    end
    return total_w
end

end
