module Algorithms

using ..CoreTypes
using Graphs
using SimpleWeightedGraphs
using LinearAlgebra
using Statistics
using DataStructures

export GreedySpanner, WGreedyWithSkeleton, RepairWithSkeleton, WAssignmentGreedy, FilteredGreedy, run_algorithm

# -----------------------------------------------------------------------------
# Algorithm Structs
# -----------------------------------------------------------------------------

struct GreedySpanner <: AbstractSpannerAlgorithm end

struct WGreedyWithSkeleton <: AbstractSpannerAlgorithm 
    skeleton_t::Float64
end
WGreedyWithSkeleton() = WGreedyWithSkeleton(1.25)

struct RepairWithSkeleton <: AbstractSpannerAlgorithm 
    skeleton_t::Float64
end
RepairWithSkeleton() = RepairWithSkeleton(1.25)

struct WAssignmentGreedy <: AbstractSpannerAlgorithm 
    skeleton_t::Float64
end
WAssignmentGreedy() = WAssignmentGreedy(1.25)

struct FilteredGreedy <: AbstractSpannerAlgorithm end

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

function get_all_edges(n, points)
    edges = Vector{Tuple{Int, Int, Float64}}(undef, n * (n - 1) ÷ 2)
    idx = 1
    for i in 1:n
        for j in i+1:n
            d = norm(points[i] - points[j])
            edges[idx] = (i, j, d)
            idx += 1
        end
    end
    sort!(edges, by = x -> x[3])
    return edges
end

function get_graph_distance(g::AbstractGraph, u::Int, v::Int, points::Vector{Point2D}, limit::Float64)
    p_end = points[v]
    h(x) = norm(points[x] - p_end)
    
    if h(u) > limit
        return Inf
    end
    
    open_set = PriorityQueue{Int, Float64}()
    enqueue!(open_set, u, h(u))
    
    g_scores = fill(Inf, nv(g))
    g_scores[u] = 0.0
    
    while !isempty(open_set)
        curr = dequeue!(open_set)
        
        if curr == v
            return g_scores[v]
        end
        
        d_curr = g_scores[curr]
        if d_curr + h(curr) > limit
            continue
        end
        
        for nbr in neighbors(g, curr)
            w = weights(g)[curr, nbr]
            tentative_g = d_curr + w
            
            if tentative_g < g_scores[nbr]
                g_scores[nbr] = tentative_g
                f_score = tentative_g + h(nbr)
                if f_score <= limit
                    open_set[nbr] = f_score
                end
            end
        end
    end
    return Inf
end

# -----------------------------------------------------------------------------
# Implementations
# -----------------------------------------------------------------------------

# 1. Greedy
function run_algorithm(algo::GreedySpanner, instance::SpannerInstance)
    start_time = time()
    n = length(instance.points)
    edges = get_all_edges(n, instance.points)
    
    g = SimpleWeightedGraph(n)
    for (u, v, dist) in edges
        limit = instance.t * dist
        d_g = get_graph_distance(g, u, v, instance.points, limit)
        if d_g > limit
            add_edge!(g, u, v, dist)
        end
    end
    
    runtime = time() - start_time
    return SpannerResult("GreedySpanner", g, runtime, Dict{Symbol, Any}())
end

# 2. w_greedy_with_skelaton
function run_algorithm(algo::WGreedyWithSkeleton, instance::SpannerInstance)
    start_time = time()
    n = length(instance.points)
    
    # 1. E <- Greedy(S, t=skeleton_t)
    edges = get_all_edges(n, instance.points) 
    g = SimpleWeightedGraph(n)
    skeleton_t = algo.skeleton_t
    for (u, v, dist) in edges
        limit = skeleton_t * dist
        d_g = get_graph_distance(g, u, v, instance.points, limit)
        if d_g > limit
            add_edge!(g, u, v, dist)
        end
    end
    
    # 2. Sort all pairs (p,q) in non-decreasing order (already in `edges`)
    # 3. For each (p,q): If w(p,q) * delta_E(p,q) > t * |pq| -> E <- E U {(p,q)}
    
    for (u, v, dist) in edges
        w_uv = instance.w_func(u, v)
        limit = (instance.t * dist) / w_uv
        
        d_g = get_graph_distance(g, u, v, instance.points, limit)
        
        if d_g > limit
            add_edge!(g, u, v, dist)
        end
    end
    
    runtime = time() - start_time
    return SpannerResult("WGreedyWithSkeleton(t_skel=$(algo.skeleton_t))", g, runtime, Dict{Symbol, Any}())
end


# 3. repair_with_skelaton
function run_algorithm(algo::RepairWithSkeleton, instance::SpannerInstance)
    start_time = time()
    n = length(instance.points)
    
    # 1. E <- Greedy(S, t=skeleton_t)
    edges = get_all_edges(n, instance.points)
    g = SimpleWeightedGraph(n)
    skeleton_t = algo.skeleton_t
    for (u, v, dist) in edges
        limit = skeleton_t * dist
        d_g = get_graph_distance(g, u, v, instance.points, limit)
        if d_g > limit
            add_edge!(g, u, v, dist)
        end
    end
    
    # 2. Add light edges: wt(e) <= wt(MST(S))/n
    mst_weight = 0.0
    ds = IntDisjointSets(n)
    edges_count = 0
    for (u, v, dist) in edges
        if !in_same_set(ds, u, v)
            union!(ds, u, v)
            mst_weight += dist
            edges_count += 1
            if edges_count == n - 1; break; end
        end
    end
    threshold = mst_weight / n
    
    for (u, v, dist) in edges
        if dist <= threshold
             add_edge!(g, u, v, dist)
        end
    end
    
    # 3. Sort points by distance to centroid
    center = mean(instance.points)
    perm = sortperm(1:n, by = i -> norm(instance.points[i] - center))
    
    # 4. While loop until valid
    max_iters = 100 
    iter = 0
    
    while true
        iter += 1
        
        added_in_this_pass = false
        
        for p_idx in 1:n
            p = perm[p_idx]
            for q_idx in 1:n
                q = perm[q_idx]
                if p == q; continue; end
                
                curr_q = q
                
                while curr_q != p
                    path_edges = a_star(g, curr_q, p) 
                    
                    if isempty(path_edges)
                        break 
                    end
                    
                    V = Vector{Int}()
                    push!(V, src(path_edges[1])) 
                    for e in path_edges
                        push!(V, dst(e))
                    end
                    
                    path_weight = 0.0
                    restarted = false
                    
                    for i in 2:length(V)
                        u_node = V[i-1]
                        v_node = V[i]
                        d_edge = weights(g)[u_node, v_node]
                        path_weight += d_edge
                        
                        v_1 = V[1]
                        v_i = V[i]
                        
                        w_val = instance.w_func(v_1, v_i)
                        d_euclid = norm(instance.points[v_1] - instance.points[v_i])
                        
                        if w_val * path_weight > instance.t * d_euclid
                            # Violation found!
                            dist_new = d_euclid 
                            add_edge!(g, v_1, v_i, dist_new)
                            
                            curr_q = v_i
                            added_in_this_pass = true
                            restarted = true
                            break 
                        end
                    end
                    
                    if !restarted
                        break 
                    end
                end
            end
        end
        
        if !added_in_this_pass
            break
        end
        if iter > max_iters
            println("RepairWithSkeleton: Max iterations reached")
            break
        end
    end
    
    runtime = time() - start_time
    return SpannerResult("RepairWithSkeleton(t_skel=$(algo.skeleton_t))", g, runtime, Dict{Symbol, Any}())
end

# 4. w_assignment_greedy
function run_algorithm(algo::WAssignmentGreedy, instance::SpannerInstance)
    start_time = time()
    n = length(instance.points)
    
    # 1. E <- Greedy(S, t=skeleton_t)
    edges = get_all_edges(n, instance.points)
    g = SimpleWeightedGraph(n)
    skeleton_t = algo.skeleton_t
    for (u, v, dist) in edges
        limit = skeleton_t * dist
        d_g = get_graph_distance(g, u, v, instance.points, limit)
        if d_g > limit
            add_edge!(g, u, v, dist)
        end
    end
    
    # 2. E' <- { (p,q) in E | w(p,q) > 1/1.5 }
    W = zeros(Float64, n, n)
    for i in 1:n, j in 1:n
        if i != j; W[i, j] = instance.w_func(i, j); end
    end
    
    e_prime = Tuple{Int, Int, Float64}[]
    for e in Graphs.edges(g)
        u, v = src(e), dst(e)
        if W[u, v] > (1.0 / skeleton_t)
             push!(e_prime, (u, v, weights(g)[u,v]))
        end
    end
    sort!(e_prime, by = x -> x[3])
    
    # 4. Propagate
    t_val = instance.t
    assignments_count = 0
    
    for (p, q, _) in e_prime
        w_pq = W[p, q]
        term = clamp(w_pq / (t_val * sqrt(2)), -1.0, 1.0)
        angle_limit = (π / 4) - asin(term)
        
        for r in neighbors(g, p)
            vec_pq = instance.points[q] - instance.points[p]
            vec_pr = instance.points[r] - instance.points[p]
            num = dot(vec_pq, vec_pr)
            den = norm(vec_pq) * norm(vec_pr)
            ang = acos(clamp(num / den, -1.0, 1.0))
            
            if ang < angle_limit
                 new_w = max(W[r, q], w_pq)
                 if new_w > W[r, q]
                     assignments_count += 1
                     W[r, q] = new_w
                     W[q, r] = new_w
                 end
            end
        end
        
        for r in neighbors(g, q)
            vec_pq = instance.points[q] - instance.points[p]
            vec_rq = instance.points[q] - instance.points[r]
            num = dot(vec_pq, vec_rq)
            den = norm(vec_pq) * norm(vec_rq)
            ang = acos(clamp(num / den, -1.0, 1.0))
            
            if ang < angle_limit
                new_w = max(W[p, r], w_pq)
                if new_w > W[p, r]
                    assignments_count += 1
                    W[p, r] = new_w
                    W[r, p] = new_w
                end
            end
        end
    end
    
    # 5. Return Weighted_Greedy_With_Skeleton(S, t, w_new)
    g_final = deepcopy(g) 
    edges = get_all_edges(n, instance.points) 
    
    for (u, v, dist) in edges
        w_uv = W[u, v]
        limit = (instance.t * dist) / w_uv
        d_g = get_graph_distance(g_final, u, v, instance.points, limit)
        if d_g > limit
            add_edge!(g_final, u, v, dist)
        end
    end
    
    runtime = time() - start_time
    stats = Dict{Symbol, Any}(:assignments_count => assignments_count)
    return SpannerResult("WAssignmentGreedy(t_skel=$(algo.skeleton_t))", g_final, runtime, stats)
end

# 5. filtered_greedy
function run_algorithm(algo::FilteredGreedy, instance::SpannerInstance)
    start_time = time()
    n = length(instance.points)
    
    # 1. E <- Greedy(S, t)
    edges = get_all_edges(n, instance.points)
    g = SimpleWeightedGraph(n)
    for (u, v, dist) in edges
        limit = instance.t * dist
        d_g = get_graph_distance(g, u, v, instance.points, limit)
        if d_g > limit
            add_edge!(g, u, v, dist)
        end
    end
    
    # Sort E from longest to shortest
    current_edges = []
    for e in Graphs.edges(g)
        u, v = src(e), dst(e)
        d = weights(g)[u, v]
        push!(current_edges, (u, v, d))
    end
    sort!(current_edges, by = x -> x[3], rev=true)
    
    # Loop
    for (u, v, dist_uv) in current_edges
        rem_edge!(g, u, v)
        
        # Check 1: w(p,q) * delta(p,q) > t * |pq|
        w_pq = instance.w_func(u, v)
        limit_pq = (instance.t * dist_uv) / w_pq
        d_g = get_graph_distance(g, u, v, instance.points, limit_pq)
        
        if d_g > limit_pq
            add_edge!(g, u, v, dist_uv)
        else
            # Else Check 2: For r in S, For s in S: if violation -> break
            is_valid = true
            
            valid_flag = Threads.Atomic{Bool}(true)
            
            Threads.@threads for r in 1:n
                if !valid_flag[]; continue; end
                
                # Dijkstra from r
                dists = dijkstra_shortest_paths(g, r).dists
                
                for s in 1:n
                    if r == s; continue; end
                    d_rs = dists[s]
                    
                    if d_rs == Inf
                        Threads.atomic_xchg!(valid_flag, false)
                        break
                    end
                    
                    d_euclid = norm(instance.points[r] - instance.points[s])
                    if d_euclid > 1e-10
                         w_rs = instance.w_func(r, s)
                         limit = (instance.t * d_euclid) / w_rs
                         if d_rs > limit + 1e-9
                             Threads.atomic_xchg!(valid_flag, false)
                             break
                         end
                    end
                end
            end
            
            if !valid_flag[]
                add_edge!(g, u, v, dist_uv)
            end
        end
    end
    
    runtime = time() - start_time
    return SpannerResult("FilteredGreedy", g, runtime, Dict{Symbol, Any}())
end

end
