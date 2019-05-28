#=NOTE:The purpose of these specifications is to create a container for holding
the building blocks of regression specifications=#


#holds a collection of specs for FM models in vector form
#NOTE: if this construct works well I may add it to Finometrics
struct FMSpecs{T<:Any}
  N::Ref{Int} #number of specs, reference construct allows for immutability
  captureresult::Function
  specnames::Vector{String}
  yspecs::Vector{Symbol}
  xspecs::Vector{FMExpr}
  xnames::Vector{Vector{Symbol}}
  withinspecs::Vector{FMExpr}
  clusterspecs::Vector{FMExpr}
  results::Vector{T}
end

#creates the bare constructor
function FMSpecs(sizehint::NInt = nothing;
    T::Type = FMLM,
    captureresult::Function = (f(m::FMLM)::T=m) #default result function is the full regression model
  )
  local specnames::Vector{String} = Vector{String}()
  local yspecs::Vector{Symbol} = Vector{Symbol}()
  local xspecs::Vector{FMExpr} = Vector{FMExpr}()
  local xnames::Vector{Vector{Symbol}} = Vector{Vector{Symbol}}()
  local withinspecs::Vector{FMExpr} = Vector{FMExpr}()
  local clusterspecs::Vector{FMExpr} = Vector{FMExpr}()
  local results::Vector{T} = Vector{T}()

  #makes minimal difference, but this could be handy
  specvector::Vector = [specnames, yspecs, xspecs, xnames, withinspecs, clusterspecs, results]
  if !isnothing(sizehint)
    (v::Vector->sizehint!(v, sizehint)).(specvector)
  end

  return FMSpecs{T}(Ref{Int}(0), captureresult,  specvector...)
end


#this applies the results function to each specification
function computeFMLMresults!(df::AbstractDataFrame, specs::FMSpecs)::Nothing

  #make sure the dimensions are corred
  @assert (specs.N[] == length(specs.specnames) &&
    specs.N[] == length(specs.yspecs) &&
    specs.N[] == length(specs.xspecs) &&
    specs.N[] == length(specs.xnames) &&
    specs.N[] == length(specs.withinspecs) &&
    specs.N[] == length(specs.clusterspecs))

  #runs the regressions
  #could conceivably parallelize this at some point
  for i::Int ∈ 1:specs.N[]
    m::FMLM = FMLM(df, specs.xspecs[i],  specs.yspecs[i],
      withinSym = specs.withinspecs[i], clusterSym = specs.clusterspecs[i],
      XNames=specs.xnames[i], YName = specs.yspecs[i])

    push!(specs.results, specs.captureresult(m))
  end

  return nothing
end

#provides a keyword access method for creating specs
function Base.push!(specs; specname::String = "($(specs.N[]))",
  yspec::Symbol = len > 1 ? specs.yspecs[end] : :Y,
  xspec::FMExpr = len > 1 ? specs.xspecs[end] : :X,
  xnames::Vector{Symbol} = len > 1 ? specs.xspecs[end] : [:intercept, :X],
  withinspec::FMExpr = nothing,
  clusterspec::FMExpr = nothing)

  #generate the specs
  push!(specs.specnames, specname)
  push!(specs.yspecs, yspec)
  push!(specs.xspecs, xspec)
  push!(specs.xnames, xnames)
  push!(specs.withinspecs, withinspec)
  push!(specs.clusterspecs, clusterspec)

  specs.N[] += 1

  return specs
end
