# multi-point crossover operator
struct MultiPointCrossover <: CrossoverOperator{2,2}
    points::UInt64

    MultiPointCrossover(numpoints::Integer=2) = new(numpoints)
end
numpoints(xover::MultiPointCrossover) = xover.points

function apply!(xover::MultiPointCrossover,
                targets::AbstractVector{<:AbstractVector},
                target_indices::AbstractVector{<:Integer},
                pop,
                parent_indices::AbstractVector{<:Integer})

    p1, p2 = map(i -> view(pop.individuals, :, i), parent_indices)
    c1, c2 = map(t -> view(t, :), targets)

    @inbounds copyto!(c1, p1)
    @inbounds copyto!(c2, p2)

    # sample xover positions uniformly at random and swap
    positions= sample(1:numdims(pop), numpoints(xover), replace=false)
    @inbounds c1[positions] = p2[positions]
    @inbounds c2[positions] = p1[positions]

    return targets
end

# mutation operator
struct BinaryFlipMutation <: MutationOperator
    ρ::Float64 # mutation rate
    BinaryFlipMutation(ρ::Float64=.01) = new(ρ)
end
ρ(mo::BinaryFlipMutation) = mo.ρ

function apply!(mo::BinaryFlipMutation,
                            target::AbstractVector{<:Integer},
                            target_index::Integer)

    amount_positions = ceil(Int, ρ(mo) * length(target))
    positions = sample(1:length(target), amount_positions, replace=false)
    @inbounds foreach(p -> target[p] = target[p] > 0 ? 0 : 1, positions)

    return target
end
