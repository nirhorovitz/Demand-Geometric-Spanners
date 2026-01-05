using Pkg
Pkg.activate(".")

using SpannerComparison
using Dates
using Printf
using Statistics
using JLD2
using Plots

function main()
    println("Initializing Parameter Tuning Experiments (Fine-grained)...")
    
    # Configuration
    N = 100 
    T_target = 1.1 
    
    # New parameter sweep values
    t_skeleton_values = [1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5]
    
    # Generate Instance
    println("Generating random instance with $N points...")
    instance = generate_random_instance(N, T_target; seed=123)
    
    results = []
    
    println("Running experiments...")
    
    for t_skel in t_skeleton_values
        println("\n--- Testing with Skeleton t = $t_skel ---")
        
        # We run WAssignmentGreedy as requested, and others for context/comparison
        algos = [
            Algorithms.WGreedyWithSkeleton(t_skel),
            Algorithms.RepairWithSkeleton(t_skel),
            Algorithms.WAssignmentGreedy(t_skel)
        ]
        
        for algo in algos
            res = Algorithms.run_algorithm(algo, instance)
            res = Analysis.compute_stats(instance, res)
            
            assign_count = get(res.stats, :assignments_count, NaN)
            
            push!(results, Dict(
                :algorithm => res.algorithm_name,
                :skeleton_t => t_skel,
                :time => res.runtime_seconds,
                :edges => res.stats[:num_edges],
                :weight_ratio => res.stats[:weight_ratio],
                :max_stretch => res.stats[:max_stretch_found],
                :valid => res.stats[:is_valid_spanner],
                :assignments => assign_count,
                # Clean name for grouping
                :algo_type => split(res.algorithm_name, "(")[1]
            ))
            
            println("  Finished $(res.algorithm_name): Time=$(round(res.runtime_seconds, digits=4))s, Edges=$(res.stats[:num_edges]), Valid=$(res.stats[:is_valid_spanner])")
        end
    end
    
    # Export JLD2
    output_dir = "experiments/parameter_tuning"
    mkpath(output_dir)
    jld_file = joinpath(output_dir, "fine_grained_sweep.jld2")
    println("\nSaving results to $jld_file...")
    save(jld_file, "results", results, "instance", instance)
    
    # --- Plotting ---
    println("Generating plots...")
    
    # Group data by algorithm
    alg_names = unique([r[:algo_type] for r in results])
    
    # 1. Edges vs t
    p1 = plot(title="Number of Edges vs Skeleton t", xlabel="Skeleton t", ylabel="Edges", legend=:topright)
    for name in alg_names
        data = filter(r -> r[:algo_type] == name, results)
        if isempty(data); continue; end
        # Sort by t just in case
        sort!(data, by = x -> x[:skeleton_t])
        
        X = [d[:skeleton_t] for d in data]
        Y = [d[:edges] for d in data]
        plot!(p1, X, Y, label=name, marker=:circle, linewidth=2)
    end
    
    # 2. Weight Ratio vs t
    p2 = plot(title="Weight Ratio vs Skeleton t", xlabel="Skeleton t", ylabel="Weight Ratio", legend=:topright)
    for name in alg_names
        data = filter(r -> r[:algo_type] == name, results)
        sort!(data, by = x -> x[:skeleton_t])
        X = [d[:skeleton_t] for d in data]
        Y = [d[:weight_ratio] for d in data]
        plot!(p2, X, Y, label=name, marker=:circle, linewidth=2)
    end
    
    # 3. Assignments vs t (Only for WAssignmentGreedy)
    p3 = plot(title="Assignments Count vs Skeleton t", xlabel="Skeleton t", ylabel="Assignments", legend=:topright)
    data_assign = filter(r -> r[:algo_type] == "WAssignmentGreedy", results)
    sort!(data_assign, by = x -> x[:skeleton_t])
    X_a = [d[:skeleton_t] for d in data_assign]
    Y_a = [d[:assignments] for d in data_assign]
    plot!(p3, X_a, Y_a, label="WAssignmentGreedy", marker=:circle, linewidth=2, color=:green)
    
    # Combine or save separate? Let's save a combined summary and individual ones.
    p_combined = plot(p1, p2, p3, layout=(3,1), size=(800, 900))
    
    png_path = joinpath(output_dir, "sweep_summary.png")
    savefig(p_combined, png_path)
    println("Saved plot to $png_path")
    
    # Print Summary Table
    println("\n" * "="^100)
    @printf("%-25s | %-6s | %-8s | %-6s | %-8s | %-8s | %-6s | %-6s\n", 
        "Algorithm", "t_skel", "Time(s)", "Edges", "WtRatio", "Stretch", "Valid", "Assign")
    println("-"^100)
    
    for r in results
        name_short = r[:algo_type]
        @printf("%-25s | %-6.2f | %-8.4f | %-6d | %-8.2f | %-8.2f | %-6s | %-6s\n", 
            name_short,
            r[:skeleton_t],
            r[:time],
            r[:edges],
            r[:weight_ratio],
            r[:max_stretch],
            string(r[:valid]),
            isnan(r[:assignments]) ? "-" : string(Int(r[:assignments]))
        )
    end
    println("="^100)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
