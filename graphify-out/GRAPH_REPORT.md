# Graph Report - .  (2026-04-08)

## Corpus Check
- Large corpus: 178 files · ~643,833 words. Semantic extraction will be expensive (many Claude tokens). Consider running on a subfolder, or use --no-semantic to run AST-only.

## Summary
- 400 nodes · 489 edges · 49 communities detected
- Extraction: 70% EXTRACTED · 30% INFERRED · 0% AMBIGUOUS · INFERRED: 148 edges (avg confidence: 0.53)
- Token cost: 0 input · 0 output

## God Nodes (most connected - your core abstractions)
1. `main()` - 13 edges
2. `dgf_one_pass()` - 11 edges
3. `RunResult` - 10 edges
4. `AlgorithmProtocol` - 10 edges
5. `_bisect()` - 8 edges
6. `_greedy_local()` - 7 edges
7. `_plot_table()` - 7 edges
8. `compute_metrics()` - 7 edges
9. `_exists_violation_pairs()` - 7 edges
10. `faster_dgf_one_pass()` - 6 edges

## Surprising Connections (you probably didn't know these)
- `Algorithm registry: dispatch by id, clear errors on unknown.` --uses--> `RunResult`  [INFERRED]
  algorithms/registry.py → core/run_result.py
- `Decorator to register an algorithm by id.` --uses--> `RunResult`  [INFERRED]
  algorithms/registry.py → core/run_result.py
- `Get algorithm by id.      Raises:         KeyError: With clear message if algori` --uses--> `RunResult`  [INFERRED]
  algorithms/registry.py → core/run_result.py
- `Run algorithm by id with full protocol signature.      Dispatches to the registe` --uses--> `RunResult`  [INFERRED]
  algorithms/registry.py → core/run_result.py
- `Return algorithms in canonical order; unknown IDs appended sorted.` --uses--> `RunResult`  [INFERRED]
  algorithms/registry.py → core/run_result.py

## Hyperedges (group relationships)
- **DGF Core Properties Triad** — readme_edge_demand_spanner, readme_sparsity, readme_lightness [EXTRACTED 1.00]
- **Algorithm Evaluation Framework** — baseline_filter_algorithm, baseline_repair_algorithm, baseline_greedy_algorithm [EXTRACTED 1.00]
- **Algorithm Comparison Study** — algo_greedy, algo_opt, algo_dgf_greedy_dgf, algo_yao_dgf, algo_yao_greedy_2 [EXTRACTED 1.00]
- **Performance Evaluation Metrics** — metric_runtime, metric_edges, metric_max_degree, metric_weight_mst [EXTRACTED 1.00]

## Communities

### Community 0 - "Experiment Runner & Config"
Cohesion: 0.04
Nodes (81): _adj_to_scipy(), algo_dgf(), algo_greedy(), algo_sqrt_greedy_dgf(), algo_yao(), algo_yao_dgf(), algo_yao_greedy_t(), _bootstrap() (+73 more)

### Community 1 - "Algorithm Protocol & Base"
Cohesion: 0.13
Nodes (20): AlgorithmProtocol, Protocol for spanner algorithms., Run the algorithm on a point set.          Args:             points: Point coord, Protocol, get(), list_algorithms(), Algorithm registry: dispatch by id, clear errors on unknown., Run algorithm by id and return RunResult with edges, metrics, max_stretch, meta. (+12 more)

### Community 2 - "Graph Metrics & APSP"
Cohesion: 0.13
Nodes (21): _absolute_weight(), apsp_add_edge(), apsp_shortest_path_after_removal(), _build_adjacency(), _compute_degrees(), _compute_max_stretch(), compute_metrics(), _euclidean_distances() (+13 more)

### Community 3 - "Filter Algorithm"
Cohesion: 0.13
Nodes (21): all_pairs_for_validation_check(), _all_pairs_shortest_paths(), _exists_violation_chunk(), _exists_violation_pairs(), _exists_violation_sequential(), _filter_fallback_reinsert_predicate(), _filter_pass(), Filter algorithm: greedy seed + remove/reinsert validation.  Intent: Remove/rein (+13 more)

### Community 4 - "DGF Acceleration (Rust/GPU)"
Cohesion: 0.16
Nodes (19): _any_violation(), _bisect(), _build_csr(), _csr_find(), _csr_remove(), _csr_restore(), faster_dgf_one_pass(), _get_n_workers() (+11 more)

### Community 5 - "Core Types & Config"
Cohesion: 0.18
Nodes (10): ExperimentConfig, PointSet, Core protocol dataclasses and types for Spanner experiments., Point set matching the protocol: set_id, points, n, source_type, generator_meta., Experiment configuration: problem_type, t, n, algorithms, seed., Validation for point sets and weight matrices., Validate point set: n consistency, points in [0,1]x[0,1], protocol compliance., Validate weight matrix for tw problem.     - shape n x n     - symmetric     - v (+2 more)

### Community 6 - "Comparison Report Generator"
Cohesion: 0.21
Nodes (11): build_comparison_report(), _compute_deltas(), _extract_run_metrics(), filter_trio_only(), SP-59: Before/after benchmark comparison report for filter, repair, weighted_gre, Filter per_scale runs to only include filter, repair, greedy (trio).      Args:, Write comparison_report.json and comparison_report.md to out_dir.      Args:, Extract runtime_ms_median, max_stretch, status, errors from run_data. (+3 more)

### Community 7 - "JSON Serialization Utils"
Cohesion: 0.2
Nodes (11): _json_default(), JSON serialization utilities for experiment outputs., Convert non-JSON-serializable objects for json.dumps., Recursively convert to JSON-serializable form (handles float inf/nan, numpy)., Convert metadata to JSON-serializable form.      Replaces non-serializable value, Serialize to JSON with numpy and common types supported.      Uses _json_default, Write JSON to file with numpy support., safe_dump() (+3 more)

### Community 8 - "Results Writer & Environment"
Cohesion: 0.24
Nodes (9): _gather_runtime_env(), Results writer for experiment outputs with schema validation., Convert RunResult (or similar) to dict for write_experiment_results.      Handle, Gather runtime environment metadata (JSON-serializable)., Convert a single algorithm run (success or failure) to JSON-serializable dict., Write manifest.json and results.json to a canonical output directory.      Args:, run_result_to_dict(), _run_to_dict() (+1 more)

### Community 9 - "Timing & Profiling"
Cohesion: 0.22
Nodes (8): measure_fn(), perf_counter(), Runtime instrumentation for algorithm and phase timing., Result of a timed operation with optional phase breakdown.      Attributes:, Merge timing metadata for inclusion in run metrics., Run a callable and return (result, runtime_ms).      Args:         fn: Callable, Return current time from perf_counter for manual timing.      Usage:         t0, TimingResult

### Community 10 - "Schema Validation"
Cohesion: 0.29
Nodes (9): _get_manifest_schema(), _get_results_schema(), _load_schema(), Schema validation hook for manifest and results before writing., Load JSON schema from schemas/ directory., Validate manifest payload against manifest schema.      Raises:         jsonsche, Validate results payload against results schema.      Raises:         jsonschema, validate_manifest() (+1 more)

### Community 11 - "Double Greedy Algorithm"
Cohesion: 0.27
Nodes (9): _all_pairs_shortest_paths(), Double greedy: two-stage approach — greedy skeleton plus weighted refinement., All-pairs shortest paths in graph with given edge lengths. Returns n x n matrix., Sort by distance ascending; ties broken by (i, j) for determinism., Sort by m(p,q) ascending; ties broken by (i, j) for determinism., Two-stage double-greedy spanner algorithm.      Stage 1 (skeleton): Build greedy, run(), _sort_by_distance_then_index() (+1 more)

### Community 12 - "Bucket Greedy Algorithm"
Cohesion: 0.27
Nodes (9): _all_pairs_shortest_paths(), _bucket_index(), Bucket greedy: partition candidates by distance buckets (alpha>1), sort within b, All-pairs shortest paths in graph with given edge lengths. Returns n x n matrix., Resolve alpha from config with safe default. alpha must be > 1., Compute bucket index for distance. Bucket k = [alpha^k, alpha^(k+1))., Bucket greedy spanner algorithm.      - Partition candidate pairs into distance, _resolve_alpha() (+1 more)

### Community 13 - "Base Algorithm Utilities"
Cohesion: 0.2
Nodes (9): progress_iter(), Shared algorithm protocol and base utilities., Wrap iterable with tqdm based on protocol config.      Defaults to showing progr, Resolve weight matrix. weight=None => ones convention.      Returns:         n x, Resolve candidate edge set for algorithm.      E_input=None => full graph candid, Resolve apsp_mode from config for incremental vs full Floyd-Warshall.      confi, resolve_apsp_mode(), resolve_candidates() (+1 more)

### Community 14 - "Assignment Greedy Algorithm"
Cohesion: 0.29
Nodes (9): _all_pairs_shortest_paths(), _greedy_skeleton(), _propagate_weights(), Weighted assignment greedy: propagation + refinement.  Intent: Propagation — spr, Weighted assignment greedy: propagation + refinement.      Stage 1 (skeleton): G, All-pairs shortest paths in graph with given edge lengths. Returns n x n matrix., Propagate weights along skeleton edges. Symmetric and bounded.      For each ske, Build greedy t-spanner skeleton ordered by Euclidean distance.     Add edge (p,q (+1 more)

### Community 15 - "Plane Spanner (Delaunay)"
Cohesion: 0.36
Nodes (8): build_graph_from_delaunay(), euclidean_dist(), find_max_stretch_pair(), generate_squashed_fibonacci_boundary(), main(), Builds a NetworkX graph from the Delaunay triangulation of the points., Finds the pair of nodes (u, v) with the maximum stretch.     Stretch = (shortest, Generates points based on a Fibonacci spiral (Golden Ratio),     squashes them i

### Community 16 - "Graph Candidate Edges"
Cohesion: 0.32
Nodes (7): full_graph_candidates(), get_candidate_edges(), Graph utilities for algorithm chaining and edge set handling., Validate and normalize edge set for n nodes.      - Raises ValueError if any nod, Resolve candidate edge set for algorithm input.      - E_input=None => full grap, Return all undirected edges as candidate set for n nodes.      Format: array of, validate_edge_set()

### Community 17 - "MST Computation"
Cohesion: 0.32
Nodes (7): _euclidean_distances(), mst_weight(), _prim_mst(), MST baseline weight for point sets., Compute pairwise Euclidean distances for points (n, 2). Returns n x n matrix., Compute MST total weight via Prim's algorithm.     dist_matrix is n x n symmetri, Compute MST weight for point set (full graph on nodes).      Args:         point

### Community 18 - "Weighted Ordered Greedy"
Cohesion: 0.32
Nodes (7): _all_pairs_shortest_paths(), Weighted ordered greedy: sort by m(p,q)=|pq|/w(p,q) non-decreasing, add edge iff, All-pairs shortest paths in graph with given edge lengths. Returns n x n matrix., Add edge iff delta_E(p,q) > t * |pq|. When delta_E=inf (disconnected), add., Weighted ordered greedy spanner algorithm.      - Process pairs sorted by m(p,q), run(), _should_add_weighted_ordered()

### Community 19 - "Repair Algorithm"
Cohesion: 0.32
Nodes (7): _greedy_skeleton(), Repair-style iterative augmentation from greedy seed.  Intent: Iterative repair, Repair-style iterative augmentation from greedy seed.      Stage 1 (seed): Greed, Build greedy t-spanner skeleton ordered by Euclidean distance.     Add edge (p,q, One repair pass: consider remaining candidates ordered by m(p,q)=|pq|/w(p,q),, _repair_pass(), run()

### Community 20 - "Two-Stage Algorithm"
Cohesion: 0.32
Nodes (7): compute_stage1_params(), _is_valid_stage1_output(), Reusable two-stage composition helper for chaining algorithms., Compute canonical stage1 parameters from (t, config) for reuse decisions.      R, Check if array is valid precomputed stage1 output: (m, 2) integer-like., Compose two algorithms: stage-1 output edges are passed to stage-2 as E_input., run_two_stage()

### Community 21 - "DGF Problem Definition"
Cohesion: 0.29
Nodes (7): Edge-Demand Function w, Demand-Geometric-Spanners Problem, Edge-Demand t-Spanner Property, Lightness O(wt(MST)), Sparsity Constraint O(n), Stretch Parameter t, SpannerComparison Julia Project

### Community 22 - "Baseline Algorithm Profiles"
Cohesion: 0.48
Nodes (7): Edge Count Metric, Filter Algorithm, Greedy Algorithm, Repair Algorithm, Stretch Verification (max_stretch), Full APSP Mode, Incremental APSP Mode

### Community 23 - "Output Path Management"
Cohesion: 0.4
Nodes (5): build_results_dir_name(), Canonical output path builder for experiment results., Build canonical output directory name.      Format: problem-{t|tw}__{set_label}_, Resolve output directory path, handling collisions.      If the canonical path a, resolve_output_path()

### Community 24 - "Greedy Spanner Algorithm"
Cohesion: 0.4
Nodes (5): Weighted greedy spanner: sort by Euclidean distance, add edge iff w(p,q)*delta_E, Add edge iff w(p,q) * delta_E(p,q) > t * |pq|. When delta_E=inf (disconnected),, Weighted greedy spanner algorithm.      - Process candidate pairs sorted by Eucl, run(), _should_add_weighted_greedy()

### Community 25 - "Yao/DGF Algorithm Variants"
Cohesion: 0.33
Nodes (6): DGF Greedy (dgf_greedy_dgf), Greedy Algorithm, Opt Algorithm, Yao DGF (yao_dgf), Yao Greedy 2 (yao_greedy_2), Runtime (seconds)

### Community 26 - "Baseline Report Writer"
Cohesion: 0.5
Nodes (3): Baseline profile artifact writer for SP-48 benchmark sweep., Write baseline_profile.json and baseline_profile.md to out_dir.      Args:, write_baseline_artifact()

### Community 27 - "Repaired Greedy Algorithm"
Cohesion: 0.5
Nodes (3): Repaired greedy: greedy stage1 (defaults) -> repair stage2 (original t/weight)., Repaired greedy: greedy(defaults) -> repair(original t/weight).      Stage 1: gr, run()

### Community 28 - "Filtered Greedy Algorithm"
Cohesion: 0.5
Nodes (3): Filtered greedy: greedy stage1 (defaults) -> filter stage2 (original t/weight)., Filtered greedy: greedy(defaults) -> filter(original t/weight).      Stage 1: gr, run()

### Community 29 - "Parameter Tuning Experiments"
Cohesion: 0.5
Nodes (4): Parameter Tuning Comparison 1.1-1.5, Graph Size (n parameter), Stretch Factor (t parameter), Parameter Sweep Summary

### Community 30 - "Algorithm Concept Comparison"
Cohesion: 0.5
Nodes (4): Filter Algorithm, Filtered Greedy Algorithm, Greedy Algorithm, Repaired Greedy Algorithm

### Community 31 - "Experiment Series Overview"
Cohesion: 1.0
Nodes (2): Algorithm Comparison Series, Yao Experiment Series

### Community 32 - "Rust Reconstructed Lib"
Cohesion: 1.0
Nodes (0): 

### Community 33 - "Rust Native DGF Lib"
Cohesion: 1.0
Nodes (0): 

### Community 34 - "Algorithm Registry Init"
Cohesion: 1.0
Nodes (0): 

### Community 35 - "NumPy Dependency"
Cohesion: 1.0
Nodes (1): NumPy

### Community 36 - "Matplotlib Dependency"
Cohesion: 1.0
Nodes (1): Matplotlib

### Community 37 - "TQDM Dependency"
Cohesion: 1.0
Nodes (1): TQDM Progress Bar Library

### Community 38 - "JSON Schema Dependency"
Cohesion: 1.0
Nodes (1): JSON Schema Validation

### Community 39 - "Strict Verification Policy"
Cohesion: 1.0
Nodes (1): Strict Verification Policy

### Community 40 - "Edge Count Metric"
Cohesion: 1.0
Nodes (1): Number of Edges

### Community 41 - "Weight/MST Metric"
Cohesion: 1.0
Nodes (1): Weight / MST

### Community 42 - "Max Degree Metric"
Cohesion: 1.0
Nodes (1): Maximum Degree

### Community 43 - "Yao Series n=200"
Cohesion: 1.0
Nodes (1): Yao Experiment Series n=200

### Community 44 - "Yao Series n=114"
Cohesion: 1.0
Nodes (1): Yao Experiment Series n=114

### Community 45 - "Yao Series n=105"
Cohesion: 1.0
Nodes (1): Yao Experiment Series n=105

### Community 46 - "Yao Series n=505"
Cohesion: 1.0
Nodes (1): Yao Experiment Series n=505

### Community 47 - "DGF Algorithm Concept"
Cohesion: 1.0
Nodes (1): DGF Algorithm

### Community 48 - "Spanner Algorithm Concept"
Cohesion: 1.0
Nodes (1): Spanner Algorithm

## Knowledge Gaps
- **172 isolated node(s):** `Faster DGF (Descending Greedy Filter) v7 =======================================`, `Build sorted CSR from boolean exists matrix.`, `Find data index for (row,col). Returns -1 if absent.`, `Set (u,v) and (v,u) to inf (Dijkstra ignores inf-weight edges). O(log degree).`, `Restore (u,v) and (v,u). O(log degree).` (+167 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Experiment Series Overview`** (2 nodes): `Algorithm Comparison Series`, `Yao Experiment Series`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Rust Reconstructed Lib`** (1 nodes): `reconstructed_lib.rs`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Rust Native DGF Lib`** (1 nodes): `lib.rs`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Algorithm Registry Init`** (1 nodes): `__init__.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `NumPy Dependency`** (1 nodes): `NumPy`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Matplotlib Dependency`** (1 nodes): `Matplotlib`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `TQDM Dependency`** (1 nodes): `TQDM Progress Bar Library`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `JSON Schema Dependency`** (1 nodes): `JSON Schema Validation`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Strict Verification Policy`** (1 nodes): `Strict Verification Policy`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Edge Count Metric`** (1 nodes): `Number of Edges`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Weight/MST Metric`** (1 nodes): `Weight / MST`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Max Degree Metric`** (1 nodes): `Maximum Degree`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Yao Series n=200`** (1 nodes): `Yao Experiment Series n=200`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Yao Series n=114`** (1 nodes): `Yao Experiment Series n=114`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Yao Series n=105`** (1 nodes): `Yao Experiment Series n=105`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Yao Series n=505`** (1 nodes): `Yao Experiment Series n=505`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `DGF Algorithm Concept`** (1 nodes): `DGF Algorithm`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Spanner Algorithm Concept`** (1 nodes): `Spanner Algorithm`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `AlgorithmProtocol` connect `Algorithm Protocol & Base` to `Base Algorithm Utilities`?**
  _High betweenness centrality (0.004) - this node is a cross-community bridge._
- **Are the 12 inferred relationships involving `main()` (e.g. with `algo_greedy()` and `algo_dgf()`) actually correct?**
  _`main()` has 12 INFERRED edges - model-reasoned connections that need verification._
- **Are the 9 inferred relationships involving `dgf_one_pass()` (e.g. with `_DgfProfiler` and `.add()`) actually correct?**
  _`dgf_one_pass()` has 9 INFERRED edges - model-reasoned connections that need verification._
- **Are the 7 inferred relationships involving `RunResult` (e.g. with `run_with_metrics()` and `Algorithm registry: dispatch by id, clear errors on unknown.`) actually correct?**
  _`RunResult` has 7 INFERRED edges - model-reasoned connections that need verification._
- **Are the 6 inferred relationships involving `AlgorithmProtocol` (e.g. with `Algorithm registry: dispatch by id, clear errors on unknown.` and `Decorator to register an algorithm by id.`) actually correct?**
  _`AlgorithmProtocol` has 6 INFERRED edges - model-reasoned connections that need verification._
- **Are the 6 inferred relationships involving `_bisect()` (e.g. with `_parallel_apsp()` and `_any_violation()`) actually correct?**
  _`_bisect()` has 6 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Faster DGF (Descending Greedy Filter) v7 =======================================`, `Build sorted CSR from boolean exists matrix.`, `Find data index for (row,col). Returns -1 if absent.` to the rest of the system?**
  _172 weakly-connected nodes found - possible documentation gaps or missing edges._