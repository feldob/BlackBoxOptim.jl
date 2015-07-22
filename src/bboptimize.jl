# Setup a fixed-dimensional problem
function setup_problem(problem::OptimizationProblem, parameters::Parameters = EMPTY_PARAMS)
  return problem, chain(DefaultParameters, parameters)
end

# Create a fixed-dimensional problem given
#   any-dimensional problem and a number of dimensions as a parameter
function setup_problem(family::FunctionBasedProblemFamily, parameters::Parameters = EMPTY_PARAMS)
  params = chain(DefaultParameters, parameters)

  # If an anydim problem was given the dimension param must have been specified.
  if params[:NumDimensions] == :NotSpecified
    throw(ArgumentError("You MUST specify the number of dimensions in a solution when a problem family is given"))
  end
  problem = fixed_dim_problem(family, parameters[:NumDimensions])

  return problem, params
end

# Create a fixed-dimensional problem given
#   a function and a search range + number of dimensions.
function setup_problem(func::Function, parameters::Parameters = EMPTY_PARAMS)
  params = chain(DefaultParameters, parameters)

  ss = check_and_create_search_space_from_parameters(params)

  # Now create an optimization problem with the given information. We currently reuse the type
  # from our pre-defined problems so some of the data for the constructor is dummy.
  problem = convert(FunctionBasedProblem, func, "", MinimizingFitnessScheme, ss) # FIXME v0.3 workaround

  # Create a random solution from the search space and ensure that the given function returns a Number.
  ind = rand_individual(BlackBoxOptim.search_space(problem))
  res = fitness(ind, problem)
  if !isa(res, Number)
    # FIXME when also supporting multi-objective opt it should not return a number
    throw(ArgumentError("The supplied function does NOT return a number when called with a potential solution (when called with $(ind) it returned $(res)) so we cannot optimize it!"))
  end

  return problem, params
end

function update_params!(parameters::Associative, kwargs)
  for (k,v) in kwargs
    parameters[k] = v
  end
end

function bboptimize(optctrl::OptController; kwargs...)
  if length(kwargs) > 0
    kwargs = convert_to_dict_symbol_any(kwargs)
    update_parameters!(optctrl, kwargs)
  end
  run!(optctrl)
end

function bboptimize(functionOrProblem, parameters::Associative = @compat(Dict{Any,Any}()); kwargs...)
  optctrl = bbsetup(functionOrProblem, parameters; kwargs...)
  run!(optctrl)
end

function convert_and_chain(parameters, kwargs)
  kwargs = convert_to_dict_symbol_any(kwargs)
  parameters = convert_to_dict_symbol_any(parameters)
  chain(parameters, kwargs)
end

function bbsetup(functionOrProblem, parameters::Associative = @compat(Dict{Any,Any}()); kwargs...)

  parameters = convert_and_chain(parameters, kwargs)
  problem, params = setup_problem(functionOrProblem, chain(DefaultParameters, parameters))
  check_valid(params)

  optimizer_func = ValidMethods[params[:Method]]
  optimizer = optimizer_func(problem, params)

  # Now set up a controller for this problem. This will handle
  # application of optimizer, checking for termination conditions
  # as well as collecting the statistics about the number of function evals,
  # keep an archive and top list of candidates.
  ctrl = OptController(optimizer, problem, params)

  return ctrl

end
