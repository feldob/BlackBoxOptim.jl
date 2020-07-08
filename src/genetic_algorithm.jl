using Distributions

const GA_DefaultOptions = ParamsDict(
    :CrossoverOperator => MultiPointCrossover(2),
    :MutationOperator => BinaryFlipMutation(.01),
    :EmbeddingOperator => NoEmbedding(),
    :Selector => SimpleSelector()
)

mutable struct GeneticAlgorithm{P<:Population,
                                S<:IndividualsSelector,
                                X<:CrossoverOperator,
                                M<:MutationOperator,
                                E<:EmbeddingOperator} <: EvolutionaryOptimizer
    name::String
    population::P
    so::S
    xo::X        # xover operator
    mo::M        # mutation operator
    eo::E         # embedding operator

    function GeneticAlgorithm(name, pop::P, so::S, xo::X, mo::M, eo::E) where {P,S,X,M,E}
        return new{P,S,X,M,E}(name, pop, so, xo, mo, eo)
    end
end

function ask(ga::GeneticAlgorithm)
    selected_indices = select(ga.so, population(ga), numparents(ga.xo) + numchildren(ga.xo))

    @inbounds parent_indices = selected_indices[1:numparents(ga.xo)]
    @inbounds target_indices = selected_indices[(numparents(ga.xo) + 1):end]

    targets = map( ti -> acquire_candi(population(ga), ti), target_indices)
    trials = map( t -> acquire_candi(population(ga), t), targets)
    trial_params = map( t -> view(t.params, :), trials)

    # xover
    apply!(ga.xo, trial_params, target_indices, population(ga), parent_indices)
    # mutations
    foreach(idx -> apply!(ga.mo, trial_params[idx], target_indices[idx]), eachindex(target_indices))

    return post_process(ga, targets, trials, ga.xo)
end

# implementation and integration into bbo

function ga(problem::OptimizationProblem,
                            options::Parameters = EMPTY_PARAMS,
                            name = "GA/ga")
    opts = chain(GA_DefaultOptions, options)
    pop = get(opts, :Population, population(problem, options))
    so = opts[:Selector]
    xo = opts[:CrossoverOperator]
    mo = opts[:MutationOperator]
    eo = opts[:EmbeddingOperator]

    return GeneticAlgorithm(name, pop, so, xo, mo, eo)
end
