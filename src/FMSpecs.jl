#=NOTE:The purpose of these specifications is to create a container for holding
the building blocks of regression specifications=#


#holds a collection of specs for FM models in vector form
#NOTE: if this construct works well I may add it to Finometrics
struct FMSpecs{T<:Any}
  N::Ref{Int} #number of specs, reference construct allows for immutability
  aggfunc::Function
  specnames::Vector{String}
  yspecs::Vector{Symbol}
  xspecs::Vector{FMExpr}
  xnames::Vector{Vector{String}}
  withinspecs::Vector{FMExpr}
  clusterspecs::Vector{Union{NSymbol, Vector{Symbol}}}
  results::Vector{T}
end



#creates the bare constructor
function FMSpecs(sizehint::NInt = nothing,
    ::Type{T} = FMLM;
    aggfunc::Function = (f(m::FMLM)::T=m) #default result function is the full regression model
  )  where T
  local specnames::Vector{String} = Vector{String}()
  local yspecs::Vector{Symbol} = Vector{Symbol}()
  local xspecs::Vector{FMExpr} = Vector{FMExpr}()
  local xnames::Vector{Vector{String}} = Vector{Vector{String}}()
  local withinspecs::Vector{FMExpr} = Vector{FMExpr}()
  local clusterspecs::Vector{Union{NSymbol, Vector{Symbol}}} = Vector{Union{NSymbol, Vector{Symbol}}}()
  local results::Vector{T} = Vector{T}()

  #makes minimal difference, but this could be handy
  specvector::Vector = [specnames, yspecs, xspecs, xnames, withinspecs, clusterspecs, results]
  if !isnothing(sizehint)
    (v::Vector->sizehint!(v, sizehint)).(specvector)
  end

  return FMSpecs{T}(Ref{Int}(0), aggfunc,  specvector...)
end

Base.show(io::IO, fs::FMSpecs) = print(io, "N: $(fs.N)
  specnames: $(fs.specnames)
  yspecs: $(fs.yspecs)
  xspecs: $(fs.xspecs)
  xnames: $(fs.xnames)
  withinspecs: $(fs.withinspecs)
  clusterspecs: $(fs.clusterspecs)
  results: $(fs.results)")

#this applies the results function to each specification
function computeFMLMresults!(dfs::S,
    specs::FMSpecs{T}, ::Type{M} = Matrix{Float64}, ::Type{V} = Vector{Float64}
    ; parallel::Bool=false, containsmissings::Bool=true, qrtype::Type = M
    )::Nothing where {S<:Union{Vector{<:AbstractDataFrame}, GroupedDataFrame},
    T<:Any, M<:AbstractMatrix, V<:AbstractVector}

  #make sure the dimensions are corred
  (specs.N[] == length(specs.specnames) &&
    specs.N[] == length(specs.yspecs) &&
    specs.N[] == length(specs.xspecs) &&
    specs.N[] == length(specs.xnames) &&
    specs.N[] == length(specs.withinspecs) &&
    specs.N[] == length(specs.clusterspecs)
    ) || error("Dimension mismatch error 43434 FMSpecs")

  (specs.N[] == length(dfs)) && "Input df views dimenion mismatch:
    dfs length: $(length(dfs)), specs.N[]: $(specs.N[])"

  results::Vector{T} = Vector{T}(undef, specs.N[])
  #runs the regressions
  #could conceivably parallelize this at some point
  if parallel #is there a problem here? if so I don't know what it is, but I have had some crashes
    Threads.@threads for i::Int ∈ 1:specs.N[]
      m::FMLM = FMLM(dfs[i], specs.xspecs[i],  specs.yspecs[i], M, V,
        withinsym = specs.withinspecs[i], clustersyms = specs.clusterspecs[i],
        Xnames=specs.xnames[i], Yname = "$(specs.yspecs[i])",
        containsmissings=containsmissings, qrtype=qrtype)
      results[i] = specs.aggfunc(m)
    end
    else
      for i::Int ∈ 1:specs.N[]
        m::FMLM = FMLM(dfs[i], specs.xspecs[i],  specs.yspecs[i], M, V,
          withinsym = specs.withinspecs[i], clustersyms = specs.clusterspecs[i],
          Xnames=specs.xnames[i], Yname = "$(specs.yspecs[i])",
          containsmissings=containsmissings, qrtype=qrtype)
        results[i] = specs.aggfunc(m)
      end
    end

  (r::T->push!(specs.results, r)).(results)

  return nothing
end

#convenience method for providing a single dataframe
function computeFMLMresults!(df::T, specs::FMSpecs,
  ::Type{M} = Matrix{Float64}, ::Type{V} = Vector{Float64};
  parallel::Bool=false, containsmissings::Bool = true, qrtype::Type=M) where{
    T<:AbstractDataFrame, M<:AbstractMatrix, V<:AbstractVector}

  return computeFMLMresults!((i::Int->view(df, :, :)).(1:(specs.N[])), specs,M,V,
    parallel=parallel, containsmissings=containsmissings, qrtype=qrtype)
end

#provides a keyword access method for creating specs
function Base.push!(specs; specname::String = "($(specs.N[]+1))",
  yspec::Symbol = len > 1 ? specs.yspecs[end] : :Y,
  xspec::FMExpr = len > 1 ? specs.xspecs[end] : :X,
  xnames::Vector{<:DField} = len > 1 ? specs.xspecs[end] : ["intercept", "X"],
  withinspec::FMExpr = nothing,
  clusterspec::Union{NSymbol, Vector{Symbol}} = nothing)

  #generate the specs
  push!(specs.specnames, specname)
  push!(specs.yspecs, yspec)
  push!(specs.xspecs, xspec)
  push!(specs.xnames, eltype(xnames) == String ? xnames : (string).(xnames))
  push!(specs.withinspecs, withinspec)
  push!(specs.clusterspecs, clusterspec)

  specs.N[] += 1

  return specs
end

Base.length(specs::FMSpecs) = specs.N[]
