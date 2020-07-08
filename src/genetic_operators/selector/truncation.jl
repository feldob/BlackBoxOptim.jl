"""
Truncation `IndividualsSelector`, with size.

"""
struct TruncationSelector <: IndividualsSelector
    k::Integer
end

#function select(::TruncationSelector, population::FitPopulation, n::Integer)
#    #TODO implement
#end
