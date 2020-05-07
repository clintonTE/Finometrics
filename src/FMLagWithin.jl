
#lags a variable in the time series
#assumes the order field is intended to be in ascending order
function lagwithin!(df::DataFrame,
  targets::Vector{<:DField},
  groups::Vector{<:DField},
  period::DField;
  lags::Int=1,
  laggednames::Vector{<:DField} = (lags ≠ 1 ?
    (s::DField->Symbol(:L, lags, s)).(targets) : (s::DField->Symbol(:L, s)).(targets)))::Nothing

  T::Type = eltype(df[!,period]) #the type of the period column

  issorted(df,[groups; period]) || error("df must be sorted in lagwithin")

  Ntargets::Int = length(targets)   #pre-allocate the space
  for t ∈ 1:Ntargets
    df[!,laggednames[t]] =
      Vector{Union{eltype(df[!,targets[t]]), Missing}}(undef, size(df, 1))
  end

  for subdf ∈ groupby(df, groups)
    Nsub::Int = size(subdf, 1)

    #only run these routines if we need to
    if Nsub > lags
      for i::Int ∈ (lags+1):Nsub #iterate for all values
        for t::Int ∈ 1:Ntargets
          cur::T = subdf[i, period]

          #first see if we can lag the easy way
          if subdf[i-lags, period] == cur - lags
            subdf[i,laggednames[t]] = subdf[i-lags, targets[t]]
          elseif lags ≠ 1 #now check the hard way (finding the lagged year), pointless if lags==1
            location::NInt = findfirst(isequal(cur-lags), subdf[1:(i-1),period])

            #record the lagged value if it is available
            if !isnothing(location)
              subdf[i,laggednames[t]] = subdf[location, targets[t]]
            end
          end
        end #targets for loop
      end #periods for loop
    end
  end

  return nothing
end

#helper method to handle the case of a single target and/or group
lagwithin!(df::DataFrame,
  targets::Union{DField, Vector{<:DField}},
  groups::Union{DField, Vector{<:DField}},
  period::DField;
  lags::Int=1,
  laggedname::NDField = nothing,
  laggednames::Vector{<:DField} = (lags ≠ 1 ?
    (s::DField->Symbol(:L, lags, s)).([targets;]) : (s::DField->Symbol(:L, s)).([targets;]))
    )::Nothing = (lagwithin!(df, [targets;], [groups;], period,
      lags=lags, laggednames=[something(laggedname, laggednames);]))


#creates a differenced column
function differencewithin!(df::DataFrame,
  targets::Vector{<:DField},
  groups::Vector{<:DField},
  period::DField;
  differencednames::Vector{<:DField} = (s::DField->Symbol(:D, s)).(targets),
  createlag::Bool=true,
  deletelag::Bool=true,
  laggednames::Vector{<:DField} = (deletelag ?
    (s::DField->Symbol(:L, s, :_temp)).(targets) : (s->Symbol(:L, s)).(targets))
  )::Nothing

  createlag && lagwithin!(df, targets, groups, period,
    laggednames=laggednames)
  for t ∈ 1:length(targets)
    df[!, differencednames[t]] = df[!, targets[t]] .- df[!, laggednames[t]]
    deletelag && select!(df, Not(laggednames[t]))
  end

  return nothing
end

#helper method to handle the case of a single target and/or group
differencewithin!(df::DataFrame,
  targets::Union{DField, Vector{<:DField}},
  groups::Union{DField, Vector{<:DField}},
  period::DField;
  differencedname::NDField = nothing,
  differencednames::Vector{<:DField} = (s::DField->Symbol(:D, s)).([targets;]),
  createlag::Bool=true,
  deletelag::Bool=true,
  laggedname::NDField = nothing,
  laggednames::Vector{<:DField} = (deletelag ?
    (s::DField->Symbol(:L, s, :_temp)).([targets;]) : (s::DField->Symbol(:L, s)).([targets;]))
  )::Nothing = differencewithin!(df, [targets;], [groups;], period,
    createlag=createlag, deletelag=deletelag,
    differencednames=[something(differencedname, differencednames);],
    laggednames=[something(laggedname, laggednames);])


##################Updated versions

###if these work, consider moving to finometrics
#lag within a group
function lagwithin2sorted(vals::AbstractVector{T},
  group::AbstractVector{G}, laggedgroup::Vector{<:Union{G, Missing}},
  ::Type{T}, ::Type{G}) where {T<:Any, G<:Any}

  lagged::Vector{Union{T, Missing}} = [missing; vals[1:(end-1)]]
  missingindicies::Vector{Bool} = ((lg,g)->
    ismissing(lg) || ismissing(g) || lg≠g).(laggedgroup,group)
  lagged[missingindicies] .= missing

  return lagged
end

function eliminatestaledates(lagged::AbstractVector{T}, date::AbstractVector{D},
  laggeddate::AbstractVector{Union{D,Missing}}, maxnotstale::Any,
  ::Type{T}, ::Type{D}) where {T<:Any, D<:Any}

  missingindicies::Vector{Bool} = ((ld, d)->
    ismissing(ld) || ismissing(d) || d - maxnotstale > ld).(laggeddate, date)
  lagged[missingindicies] .= missing

  return lagged
end

function lagwithin2sorted(
  vals::AbstractVector,
  group::AbstractVector;
  date::Union{AbstractVector{D},Nothing} = nothing,
  laggeddate::Union{Vector{Union{D,Missing}}, Nothing} = nothing,
  maxnotstale::Any = nothing,
  laggedgroup::Vector{Union{G, Missing}} =
    Vector{Union{eltype(group), Missing}}([missing; group[1:(end-1)]])) where {G<:Any, D<:Any}


  local lagged = lagwithin2sorted(vals, group, laggedgroup, eltype(vals), eltype(group))

  if !isnothing(maxnotstale)
    if isnothing(laggeddate)
      laggeddate = [missing; [date[1:(end-1)]]]
    end
    lagged = eliminatestaledates(lagged, date, laggeddate, maxnotstale, eltype(lagged), eltype(date))
  end

  return lagged
end

#lags multiple columns within a dataframe
function lagwithin2sorted!(df::DataFrame, vals::Vector{<:DField}, group::DField;
  date::NDField = nothing,
  maxnotstale::Any = nothing,
  laggedvals::Vector{<:DField} = (s->Symbol(:L, s)).(vals))


  local laggedgroup::Vector{Union{eltype(df[!,group]), Missing}}
  local laggeddate::(isnothing(date) ? Nothing : Vector{Union{eltype(df[!,date]), Missing}})


  #lag the groups, and dates if needed, once
  laggedgroup = Vector{Union{eltype(group), Missing}}([missing; df[!, group][1:(end-1)]])
  #println(!isnothing(date))
  if (!isnothing(date))
    laggeddate = [missing; df[!,date][1:(end-1)]]
  else
    laggeddate = nothing
  end
  #println(typeof(laggeddate))
  #allocate the taskes for multithreading
  tasks::Vector{Task} = Vector{Task}(undef, length(vals))

  #println(typeof(laggeddate))
  #deploy the lags
  tasks = @sync((s::DField -> Threads.@spawn lagwithin2sorted(
    df[!, s], df[!, group], date = isnothing(date) ? nothing : df[!,date],
    laggeddate = laggeddate, maxnotstale = maxnotstale, laggedgroup = laggedgroup)).(vals))


  for i ∈ 1:length(vals)
    #println(typeof(@fetch tasks[i]))
    df[!, laggedvals[i]] = fetch(tasks[i])
    #=df[!, laggedvals[i]] = lagwithin2sorted(
      df[!, vals[i]], df[!, group], date = isnothing(date) ? nothing : df[!,date],
      laggeddate = convert(Union{eltype(date), Missing}, laggeddate),
      maxnotstale = maxnotstale, laggedgroup = laggedgroup)=#
  end

  return nothing
end


function lagwithin2!(df::DataFrame, vals::Vector{<:DField}, group::DField;
  date::NDField = nothing,
  maxnotstale::Any = nothing,
  laggedvals::Vector{<:DField} = (s->Symbol(:L, s)).(vals))

  if !isnothing(date)
    issorted(df, [group, date]) || error("df must be sorted by $group, $date")
  else
    issorted(df, [group]) || error("df must be sorted by $group")
  end

  return lagwithin2sorted!(df, vals, group, date=date, maxnotstale=maxnotstale, laggedvals=laggedvals)
end

function differencewithin2!(df::DataFrame,
  vals::Vector{<:DField},
  group::DField;
  date::NDField = nothing,
  differencedvals::Vector{<:DField} = (s::DField->Symbol(:D, s)).(vals),
  createlag::Bool=true,
  deletelag::Bool=true,
  maxnotstale::Any = nothing,
  laggedvals::Vector{<:DField} = (deletelag ?
    (s::DField->Symbol(:L, s, :_temp)).(vals) : (s->Symbol(:L, s)).(vals))
  )::Nothing

  createlag && lagwithin2!(df, vals, group, date=date,
    laggedvals=laggedvals, maxnotstale=maxnotstale)
  for t ∈ 1:length(vals)
    df[!, differencedvals[t]] = df[!, vals[t]] .- df[!, laggedvals[t]]
    deletelag && select!(df, Not(laggedvals[t]))
  end

  return nothing
end
