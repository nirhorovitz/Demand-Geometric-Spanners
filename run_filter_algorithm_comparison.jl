using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using SpannerComparison
using SpannerComparison.CoreTypes
using SpannerComparison.Algorithms
using SpannerComparison.Generators
using SpannerComparison.Analysis
using SpannerComparison.Visualization
using SpannerComparison.IOUtils

using Dates
using Printf
using Plots

function draw_table_image(results, N, T, seed, filepath)
    headers = ["Algorithm", "Edges", "Weight", "WtRatio", "MaxDeg", "AvgDeg", "MaxStrch", "Valid", "Time(s)"]
    data = Matrix{String}(undef, length(results), length(headers))

    for (i, res) in enumerate(results)
        s = res.stats
        data[i, 1] = res.algorithm_name
        data[i, 2] = string(get(s, :num_edges, ""))
        data[i, 3] = @sprintf("%.2f", get(s, :total_weight, NaN))
        data[i, 4] = @sprintf("%.2f", get(s, :weight_ratio, NaN))
        data[i, 5] = string(get(s, :max_degree, ""))
        data[i, 6] = @sprintf("%.2f", get(s, :avg_degree, NaN))
        data[i, 7] = @sprintf("%.4f", get(s, :max_stretch_found, NaN))
        data[i, 8] = string(get(s, :is_valid_spanner, "N/A"))
        data[i, 9] = @sprintf("%.4f", res.runtime_seconds)
    end

    col_widths = [2.5, 1.2, 1.5, 1.2, 1.2, 1.2, 1.6, 1.0, 1.5]
    total_width = sum(col_widths)

    x_centers = Float64[]
    current_x = 0.0
    for w in col_widths
        push!(x_centers, current_x + w / 2)
        current_x += w
    end

    n_rows = length(results)
    f_title = font(16, :bold)
    f_head = font(12, :bold)
    f_cell = font(12)

    plot_height = 150 + n_rows * 35
    p = plot(
        axis = false,
        grid = false,
        ticks = false,
        border = :none,
        size = (1400, plot_height),
        xlim = (0, total_width),
        ylim = (-n_rows - 1.5, 1),
    )

    title_str = "N=$N | t=$T | seed=$seed"
    annotate!(p, total_width / 2, 0.5, text(title_str, f_title, :black, :center))

    for j in 1:length(headers)
        annotate!(p, x_centers[j], -0.5, text(headers[j], f_head, :black, :center))
    end

    for i in 1:n_rows
        for j in 1:length(headers)
            if j == 1
                x_pos = x_centers[j] - col_widths[j] / 2 + 0.2
                annotate!(p, x_pos, -0.5 - i, text(data[i, j], f_cell, :black, :left))
            else
                annotate!(p, x_centers[j], -0.5 - i, text(data[i, j], f_cell, :black, :center))
            end
        end
    end

    savefig(p, filepath)
end

function parse_t_values(args)
    default_ts = [1.05, 1.1, 1.2, 1.25, 1.4, 1.5, 1.75, 2.0]
    if length(args) < 2
        return default_ts
    end
    t_arg = args[2]
    if occursin(",", t_arg)
        return [parse(Float64, strip(x)) for x in split(t_arg, ",")]
    else
        return [parse(Float64, t_arg)]
    end
end

function main()
    gr()  # ensure GR backend for consistent image output

    # Defaults
    N = 300
    seed = 42
    t_values = parse_t_values(ARGS)

    if length(ARGS) >= 1
        N = parse(Int, ARGS[1])
    end
    # ARGS[2] handled by parse_t_values
    if length(ARGS) >= 3
        seed = parse(Int, ARGS[3])
    end

    println("Running filter algorithm comparison with N=$N, t set=$(t_values), seed=$seed")

    # Fix the point set and weights once for all t to ensure a fair comparison
    base_instance = generate_random_instance(N, t_values[1]; seed=seed)

    algo_sequence = Any[
        Algorithms.GreedySpanner(),     # baseline
        Algorithms.SqrtGreedyFilter(),  # Algorithm A
    ]

    base_dir = "/Users/nhorovitz/Documents/univ_stuff/spanner project/jolia/Experiments"
    root_output_dir = joinpath(base_dir, "n=$(N)_t=$(t_values[1])")
    mkpath(root_output_dir)

    for T in t_values
        println("\n=== Running t=$T ===")
        instance = SpannerInstance(base_instance.points, base_instance.w_func, T)

        results = SpannerResult[]
        for algo in algo_sequence
            println("  -> Running $(algo)...")
            res = Algorithms.run_algorithm(algo, instance)
            res = Analysis.compute_stats(instance, res)
            push!(results, res)
            println("     Completed $(res.algorithm_name): edges=$(res.stats[:num_edges]), weight=$(round(res.stats[:total_weight], digits=3)), valid=$(res.stats[:is_valid_spanner]), time=$(round(res.runtime_seconds, digits=3))s")
        end

        output_dir = joinpath(root_output_dir, "t=$(T)")
        mkpath(output_dir)

        data_path = joinpath(output_dir, "spanner_data.jld2")
        save_results(data_path, instance, results)

        table_path = joinpath(output_dir, "summary_table.png")
        draw_table_image(results, N, T, seed, table_path)

        plot_path = joinpath(output_dir, "comparison_output.png")
        comparison_plot = visualize_results(instance, results)
        savefig(comparison_plot, plot_path)

        println("Saved outputs to $output_dir")
        println("  - data: $data_path")
        println("  - table: $table_path")
        println("  - plot:  $plot_path")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
