using BlackBoxOptim

# bbo is now extended to support not only AbstractVector{Float64} represetations of solutions, but all AbstractVector{<:Real} solution types.
# To illustrate this you here find a variant with binary representations through UInt8, instead of the default (Individual = Vector{Float64}).
# The toy problem implemented below is maxone, which maximizes the amount of ones in a binary array (as a proof of concept).
# The example uses a new optimizer, a Genetic Algorithm (GA), which by default uses the newly implemented operators BinaryFlipMutation and MultiPointCrossover.

# fitness function
maxones(genotype::AbstractVector{<:Integer})::Float64 = sum(genotype)

# search space with unsigned integers
ss = RectSearchSpace(1000, (0x00, 0x01), UInt8)

bboptimize(maxones; SearchSpace = ss,
                    MaxTime = 10,
                    FitnessScheme = MaximizingFitnessScheme,
                    Method = :ga)

# Per default, the MultiPointCrossover uses 2-point crossover and a mutation rate of 1%, which has the algorithm converge slowly. Below we set the number of points to 50 and mutation rate to .1 % to speed-up convergence.
using BlackBoxOptim: MultiPointCrossover, BinaryFlipMutation

bboptimize(maxones; CrossoverOperator = MultiPointCrossover(50),
                    MutationOperator = BinaryFlipMutation(.001),
                    SearchSpace = ss,
                    MaxTime = 10,
                    FitnessScheme = MaximizingFitnessScheme,
                    Method = :ga)
