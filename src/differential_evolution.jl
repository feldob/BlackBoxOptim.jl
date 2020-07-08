using Distributions

const DE_DefaultOptions = chain(DEX_DefaultOptions, ParamsDict(
    :SamplerRadius => 8,
))

mutable struct DiffEvoOpt{P<:Population,
                          S<:IndividualsSelector,
                          X<:GeneticOperator,
                          E<:EmbeddingOperator} <: EvolutionaryOptimizer
    # TODO when sampler and bound would be parameterized, name is no longer required -- as everything is seen in the type name
    name::String
    population::P
    so::S        # selection operator
    xo::X        # genetic operator
    eo::E        # embedding operator

    function DiffEvoOpt(name, pop::P, so::S, xo::X, eo::E) where {P,S,X,E}
        return new{P,S,X,E}(name, pop, so, xo, eo)
    end
end

ask(de::DiffEvoOpt) = evolve(de, de.xo)

function evolve(de::DiffEvoOpt, xo::GeneticOperatorsMixture)
    sel_xo, tag = next(xo)
    candidates = evolve(de, sel_xo) # use the selected operator of the mixture
    # override what was set by sel_op so that adjust!() gets called for the mixture op
    for candi in candidates
        candi.extra = xo
        candi.tag = tag
    end
    return candidates
end

function evolve(de::DiffEvoOpt, mutate::MutationOperator)
    target_index = select(de.so, population(de), 1)[1]
    target = acquire_candi(population(de), target_index)
    trial = acquire_candi(population(de), target)
    apply!(mutate, trial.params, target_index)
    return post_process(de, target, trial, mutate)
end

function evolve(de::DiffEvoOpt, crossover::CrossoverOperator)
    # Sample parents and target
    indices = select(de.so, population(de), 1 + numparents(crossover))
    parent_indices = indices[1:numparents(crossover)]
    #println("parent_indices = $(parent_indices)")
    target_index = indices[end]
    target = acquire_candi(population(de), target_index)
    trial = acquire_candi(population(de), target)
    #println("target[$(target_index)] = $(target)")

    # Crossover parents and target
    @assert numchildren(crossover)==1
    apply!(crossover, trial.params, target_index,
           population(de), parent_indices)
    return post_process(de, target, trial, crossover)
end

# Now we can create specific DE optimizers that are commonly used in the
# literature.

function diffevo(problem::OptimizationProblem, options::Parameters, name::String,
                 so::IndividualsSelector = SimpleSelector(),
                 co::DiffEvoCrossoverOperator = DiffEvoRandBin1(chain(DE_DefaultOptions, options)))
    opts = chain(DE_DefaultOptions, options)
    pop = population(problem, opts)
    DiffEvoOpt(name, pop, so, co,
               RandomBound(search_space(problem)))
end

"""
The most used `DE/rand/1/bin` variant of differential evolution.
"""
de_rand_1_bin(problem::OptimizationProblem,
              options::Parameters = EMPTY_PARAMS,
              name = "DE/rand/1/bin") = diffevo(problem, options, name)

de_rand_2_bin(problem::OptimizationProblem,
              options::Parameters = EMPTY_PARAMS,
              name = "DE/rand/2/bin") = diffevo(problem, options, name,
                                                SimpleSelector(),
                                                DiffEvoRandBin2(chain(DE_DefaultOptions, options)))

"""
The most used `DE/rand/1/bin` variant with "local geography" via radius-limited sampling.
"""
de_rand_1_bin_radiuslimited(problem::OptimizationProblem,
                            options::Parameters = EMPTY_PARAMS,
                            name = "DE/rand/1/bin/radiuslimited") =
    diffevo(problem, options, name,
            RadiusLimitedSelector(chain(DE_DefaultOptions, options)[:SamplerRadius]))

de_rand_2_bin_radiuslimited(problem::OptimizationProblem,
                            options::Parameters = EMPTY_PARAMS,
                            name = "DE/rand/2/bin/radiuslimited") =
    diffevo(problem, options, name,
            RadiusLimitedSelector(chain(DE_DefaultOptions, options)[:SamplerRadius]),
            DiffEvoRandBin2(chain(DE_DefaultOptions, options)))
