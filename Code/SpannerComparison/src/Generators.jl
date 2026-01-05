module Generators

using ..CoreTypes
using StaticArrays
using Random
using LinearAlgebra

export generate_random_instance

"""
    generate_random_instance(n::Int, t::Float64; seed::Union{Int, Nothing}=nothing)

Generates a random SpannerInstance with `n` points in the unit square.
`w` function is also randomized.
"""
function generate_random_instance(n::Int, t::Float64; seed::Union{Int, Nothing}=nothing)
    if !isnothing(seed)
        Random.seed!(seed)
    end
    
    points = [SVector{2, Float64}(rand(), rand()) for _ in 1:n]
    
    # Define a random w function. 
    # For efficiency and reproducibility, we can pre-generate weights or use a deterministic hash-based function.
    # Here, let's create a symmetric matrix of random weights to simulate a general function w(p,q).
    # Since w is from (0, 1], we verify that.
    
    # Option 1: Store a matrix (Memory intensive for N=8000 -> 64MB floats, acceptable).
    # This ensures O(1) access and consistency.
    w_matrix = rand(n, n) # (0, 1)
    # Symmetrize and ensure no zeros if needed (though (0,1] is usually > 0)
    # range (0, 1] means we might want to shift rand() output if strictly > 0 is required.
    # rand() is [0, 1). let's do (rand() + eps()) to be safe or just 1.0 - rand().
    for i in 1:n
        w_matrix[i, i] = 1.0 # w(p,p) usually doesn't matter or is 1
        for j in (i+1):n
            val = rand()
            val = val == 0.0 ? 0.0001 : val # Ensure strictly positive if required by (0, 1]
            w_matrix[i, j] = val
            w_matrix[j, i] = val
        end
    end
    
    w_func = MatrixWFunc(w_matrix)
    
    return SpannerInstance(points, w_func, t)
end

end

