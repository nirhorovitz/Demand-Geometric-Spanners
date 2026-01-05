# SpannerComparison

A Julia project to compare algorithms for building $(w,t)$-spanners.

## Setup

1. Install Julia (1.8+ recommended).
2. Clone this repository.
3. Instantiate the project dependencies:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

## Usage

Run the comparison script:

```bash
julia --project=. --threads auto scripts/run_comparison.jl
```

Note: Use `--threads auto` (or a specific number like `--threads 4`) to enable parallel processing for spanner verification, which significantly speeds up analysis for large N (1000-8000).

## Project Structure

- `src/`: Source code.
  - `SpannerComparison.jl`: Main module.
  - `Core.jl`: Type definitions (`SpannerInstance`, `SpannerResult`).
  - `Generators.jl`: Random instance generation.
  - `Algorithms.jl`: Spanner algorithms (currently stubs).
  - `Analysis.jl`: Statistical analysis and verification.
  - `Visualization.jl`: Plotting code.
  - `IO.jl`: Saving/Loading results.
- `scripts/`: Executable scripts.

## Adding Algorithms

Implement new algorithms in `src/Algorithms.jl` by defining a struct `<: AbstractSpannerAlgorithm` and extending `run_algorithm`.

## Output

- `comparison_output.png`: Visualization of the generated spanners.
- `spanner_data.jld2`: Saved data containing the instance and results.

