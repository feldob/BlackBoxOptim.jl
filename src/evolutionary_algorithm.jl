abstract type EvolutionaryOptimizer <: PopulationOptimizer end

# trace current optimization state,
# Called by OptRunController trace_progress()
function trace_state(io::IO,
                     eo::EvolutionaryOptimizer,
                     mode::Symbol)
    if mode == :compact
        return
    end

    println(io, "$(eo.name) modify state:")
    trace_state(io, eo.xo, mode)
end

function tell!(eo::EvolutionaryOptimizer,
               # archive::Archive, # Skip for now
               rankedCandidates::Vector{Candidate{F}}) where F
    n_acceptable_candidates = length(rankedCandidates)÷2
    num_better = 0
    for i in eachindex(rankedCandidates)
        candi = rankedCandidates[i]
        # accept the modified individuals from the top ranked half
        if i <= n_acceptable_candidates
            is_improved = candi.params != viewer(population(eo), candi.index)
            adjust!(eo, candi, is_improved)
        else
            is_improved = false
        end
        if is_improved
            num_better += 1
            #print("candidate = "); show(candidate); println("")
            #print("index = "); show(index); println("")
            #print("target = "); show(eo.population[index,:]); println("")
            #old = eo.population[:,index]
            accept_candi!(eo.population, candi)
        else
            # just return candidate to the pool
            release_candi(eo.population, candi)
        end
    end
    num_better
end

"""
Adjust the parameters of the method after the candidate evaluation.
"""
function adjust!(eo::EvolutionaryOptimizer,
                 candi::Candidate,
                 is_improved::Bool)
    # adjust the parameters of the operation
    old_fitness = fitness(population(eo), candi.index)
    if isnan(old_fitness)
        old_fitness = candi.fitness
    end

    adjust!(candi.extra::GeneticOperator, candi.tag, candi.index, candi.fitness,
            fitness(population(eo), candi.index), is_improved)
end

"""
Post-process the evolved pairs.
"""
function post_process(eo::EvolutionaryOptimizer,
                      targets::AbstractVector{Candidate{F}},
                      trials::AbstractVector{Candidate{F}},
                      go::GeneticOperator,
                      tag::Int = 0) where F
    pairs = map(idx -> post_process(eo, targets[idx], trials[idx], go, tag), eachindex(targets))
    return reduce(vcat, pairs)
end

function post_process(eo::EvolutionaryOptimizer,
                      target::Candidate{F},
                      trial::Candidate{F},
                      go::GeneticOperator,
                      tag::Int = 0) where F
    # embed the trial parameter vector into the search space
    apply!(eo.eo, trial.params, eo.population, target.index)
    target.extra = trial.extra = go
    target.tag = trial.tag = tag
    if trial.params != target.params
        reset_fitness!(trial, eo.population)
    end

    return Candidate{F}[trial, target]
end
