include("../src/BlackBoxOptim.jl")
using .BlackBoxOptim
using Test

const NumTestRepetitions = 100

random_candidate(n, lo, hi, I=Float64) =
    BlackBoxOptim.Candidate{Float64, I}([clamp(rand(I), lo, hi) for _ in 1:n])
