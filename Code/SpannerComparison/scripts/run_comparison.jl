using Pkg
# Activate project environment relative to this script
Pkg.activate(joinpath(@__DIR__, ".."))
# Pkg.instantiate() # Uncomment if dependencies are missing

using SpannerComparison
using Dates
using Printf
using Plots

function main()
    println("Initializing Spanner Comparison Project...")
    
    # 1. Configuration
    # Default values
    N = 50 
    T = 1.5
    
    # Parse command line arguments if provided
    # Usage: julia run_comparison.jl [N] [T]
    if length(ARGS) >= 1
        N = parse(Int, ARGS[1])
    end
    if length(ARGS) >= 2
        T = parse(Float64, ARGS[2])
    end
    
    println("Configuration: N=$N, t=$T")

    
    # 2. Generate Instance
    println("Generating random instance with $N points...")
    instance = generate_random_instance(N, T; seed=42)
    
    # 3. Run Algorithms
    algorithms = [
        Algorithms.GreedySpanner(),
        Algorithms.WGreedyWithSkeleton(),
        Algorithms.RepairWithSkeleton(),
        Algorithms.WAssignmentGreedy(),
        Algorithms.FilteredGreedy()
    ]
    
    results = SpannerResult[]
    
    println("Running algorithms...")
    for algo in algorithms
        # println("  Running $(typeof(algo))...")
        res = Algorithms.run_algorithm(algo, instance)
        
        # 4. Analyze
        # println("  Analyzing results for $(res.algorithm_name)...")
        res = Analysis.compute_stats(instance, res)
        push!(results, res)
        
        println("  Finished $(res.algorithm_name): Time=$(round(res.runtime_seconds, digits=4))s, Valid=$(res.stats[:is_valid_spanner])")
    end
    
    # 5. Visualize
    println("Visualizing...")
    
    # Prepare output directory
    output_dir = "n=$(N)_t=$(T)"
    mkpath(output_dir)
    println("Saving results to directory: $output_dir")
    
    try
        plt = Visualization.visualize_results(instance, results)
        plot_path = joinpath(output_dir, "comparison_output.png")
        savefig(plt, plot_path)
        println("  Saved plot to '$plot_path'")
    catch e
        println("  Visualization failed: $e")
    end
    
    # 6. Export
    println("Exporting data...")
    data_path = joinpath(output_dir, "spanner_data.jld2")
    IOUtils.save_results(data_path, instance, results)
    
    println("Done!")
    
    # Print summary table
    @printf("\n%-20s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s\n", 
        "Algorithm", "Time(s)", "Edges", "MaxDeg", "WtRatio", "MaxStrch", "Valid")
    println("-"^90)
    for res in results
        @printf("%-20s | %-8.4f | %-8d | %-8d | %-8.2f | %-8.2f | %-8s\n", 
            res.algorithm_name, 
            res.runtime_seconds, 
            res.stats[:num_edges], 
            res.stats[:max_degree], 
            get(res.stats, :weight_ratio, NaN),
            get(res.stats, :max_stretch_found, NaN),
            get(res.stats, :is_valid_spanner, "N/A")
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
