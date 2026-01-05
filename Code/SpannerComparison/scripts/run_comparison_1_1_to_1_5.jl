using Pkg
Pkg.activate(".")

using SpannerComparison
using Dates
using Printf
using Statistics
using JLD2
using Plots
using Graphs
using SimpleWeightedGraphs

function are_graphs_identical(g1, g2)
    if nv(g1) != nv(g2) || ne(g1) != ne(g2)
        return false
    end
    for e in edges(g1)
        u, v = src(e), dst(e)
        w1 = weight(e)
        # Check if edge exists in g2
        # SimpleWeightedGraphs doesn't have has_edge(g, u, v) optimized always, but it works.
        # weights(g)[u,v] is 0 if no edge.
        w2 = weights(g2)[u, v]
        if w2 == 0 || abs(w1 - w2) > 1e-9
            return false
        end
    end
    return true
end

function main()
    println("Initializing Experiment (Comparison 1.1 - 1.5)...")
    
    # Configuration
    N = 100 
    T_target = 1.1 
    
    # Range 1.1 to 1.5
    t_skeleton_values = collect(1.1:0.05:1.5)
    
    # Generate Instance
    println("Generating random instance with $N points...")
    instance = generate_random_instance(N, T_target; seed=123)
    
    results = []
    
    println("Running experiments...")
    
    for t_skel in t_skeleton_values
        println("\n--- Testing with Skeleton t = $(round(t_skel, digits=2)) ---")
        
        # Run WGreedyWithSkeleton
        algo_wgs = Algorithms.WGreedyWithSkeleton(t_skel)
        res_wgs = Algorithms.run_algorithm(algo_wgs, instance)
        res_wgs = Analysis.compute_stats(instance, res_wgs)
        
        # Run RepairWithSkeleton
        algo_rep = Algorithms.RepairWithSkeleton(t_skel)
        res_rep = Algorithms.run_algorithm(algo_rep, instance)
        res_rep = Analysis.compute_stats(instance, res_rep)
        
        # Run WAssignmentGreedy
        algo_wag = Algorithms.WAssignmentGreedy(t_skel)
        res_wag = Algorithms.run_algorithm(algo_wag, instance)
        res_wag = Analysis.compute_stats(instance, res_wag)
        
        # Check Identity
        identical = are_graphs_identical(res_wgs.graph, res_wag.graph)
        
        # Store Results
        # Helper to push
        function add_res(r, is_wgs_wag_identical=nothing)
            push!(results, Dict(
                :algorithm => r.algorithm_name,
                :short_name => split(r.algorithm_name, "(")[1],
                :skeleton_t => t_skel,
                :edges => r.stats[:num_edges],
                :weight_ratio => r.stats[:weight_ratio],
                :max_stretch => r.stats[:max_stretch_found],
                :assignments => get(r.stats, :assignments_count, NaN),
                :identical_to_wgs => is_wgs_wag_identical
            ))
        end
        
        add_res(res_wgs)
        add_res(res_rep)
        add_res(res_wag, identical)
        
        println("  WGreedyWithSkeleton: Edges=$(res_wgs.stats[:num_edges])")
        println("  RepairWithSkeleton:  Edges=$(res_rep.stats[:num_edges])")
        println("  WAssignmentGreedy:   Edges=$(res_wag.stats[:num_edges]), Assign=$(get(res_wag.stats, :assignments_count, 0))")
        println("  WGS == WAG?          $identical")
    end
    
    # Export
    output_dir = "experiments/parameter_tuning"
    mkpath(output_dir)
    
    # Use timestamp to avoid overwriting
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    jld_file = joinpath(output_dir, "comparison_1.1_1.5_$(timestamp).jld2")
    save(jld_file, "results", results, "instance", instance)
    
    # --- Plotting ---
    println("Generating plots...")
    
    # Titles with info
    common_title_suffix = " (N=$N, T=$T_target)"
    
    # 1. Edges vs t
    p1 = plot(title="Edges vs t" * common_title_suffix, xlabel="Skeleton t", ylabel="Edges", legend=:topright)
    
    # 2. Weight Ratio vs t
    p2 = plot(title="Weight Ratio vs t" * common_title_suffix, xlabel="Skeleton t", ylabel="Ratio", legend=:topright)
    
    # 3. Assignments (WAG only)
    p3 = plot(title="Assignments (WAG) vs t" * common_title_suffix, xlabel="Skeleton t", ylabel="Count", legend=:topright)
    
    # Extract series
    t_vals = t_skeleton_values
    
    function get_series(name, key)
        vals = []
        for t in t_vals
            # Find matching result (assuming one per t per algo)
            r = findfirst(x -> x[:short_name] == name && abs(x[:skeleton_t] - t) < 1e-5, results)
            if !isnothing(r)
                push!(vals, results[r][key])
            else
                push!(vals, NaN)
            end
        end
        return vals
    end
    
    # Plot WGreedyWithSkeleton
    y_edges_wgs = get_series("WGreedyWithSkeleton", :edges)
    y_wr_wgs = get_series("WGreedyWithSkeleton", :weight_ratio)
    plot!(p1, t_vals, y_edges_wgs, label="WGreedyWithSkel", marker=:circle, lw=2)
    plot!(p2, t_vals, y_wr_wgs, label="WGreedyWithSkel", marker=:circle, lw=2)
    
    # Plot RepairWithSkeleton
    y_edges_rep = get_series("RepairWithSkeleton", :edges)
    y_wr_rep = get_series("RepairWithSkeleton", :weight_ratio)
    plot!(p1, t_vals, y_edges_rep, label="RepairWithSkel", marker=:square, lw=2, linestyle=:dash)
    plot!(p2, t_vals, y_wr_rep, label="RepairWithSkel", marker=:square, lw=2, linestyle=:dash)
    
    # Plot WAssignmentGreedy
    y_edges_wag = get_series("WAssignmentGreedy", :edges)
    y_wr_wag = get_series("WAssignmentGreedy", :weight_ratio)
    y_assign = get_series("WAssignmentGreedy", :assignments)
    
    plot!(p1, t_vals, y_edges_wag, label="WAssignmentGreedy", marker=:diamond, markersize=6, lw=1.5, linestyle=:dot)
    plot!(p2, t_vals, y_wr_wag, label="WAssignmentGreedy", marker=:diamond, markersize=6, lw=1.5, linestyle=:dot)
    
    plot!(p3, t_vals, y_assign, label="Assignments", marker=:star, lw=2, color=:purple)
    
    p_combined = plot(p1, p2, p3, layout=(3,1), size=(800, 1000))
    
    png_path = joinpath(output_dir, "comparison_1.1_1.5_$(timestamp).png")
    savefig(p_combined, png_path)
    println("Saved plot to $png_path")
    
    # Table
    println("\n" * "="^110)
    @printf("%-20s | %-5s | %-6s | %-8s | %-8s | %-6s\n", 
        "Algorithm", "t", "Edges", "WtRatio", "Assign", "WGS==WAG")
    println("-"^110)
    
    for r in results
        identical_str = isnothing(r[:identical_to_wgs]) ? "-" : string(r[:identical_to_wgs])
        assign_str = isnan(r[:assignments]) ? "-" : string(Int(r[:assignments]))
        
        @printf("%-20s | %-5.2f | %-6d | %-8.2f | %-8s | %-6s\n", 
            r[:short_name],
            r[:skeleton_t],
            r[:edges],
            r[:weight_ratio],
            assign_str,
            identical_str
        )
    end
    println("="^110)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
