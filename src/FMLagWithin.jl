
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
  group::AbstractVector{G}, laggedgroup::Vector{<:Union{G, Missing}}, ::Type{T}, ::Type{G};
  ) where {T<:Any, G<:Any}


  #idea here is to lag everyting, then account for the boundaries between groups
  #lag everything
  lagged::Vector{Union{T, Missing}} = [missing; vals[1:(end-1)]]

  #fixes lags across groups
  missingindicies::Vector{Bool} = ((lg,g)->
    ismissing(lg) || ismissing(g) || lg≠g).(laggedgroup,group)
  lagged[missingindicies] .= missing

  return lagged
end

#no dataframe version - just lags the vector
#not called in the dataframe versions
function lagwithin2sorted(vals::TV,group::TG) where {TV, TG}
  laggedgroup = Vector{Union{eltype(group), Missing}}([missing; group[1:(end-1)]])
  return lagwithin2sorted(vals, group, laggedgroup, eltype(TV), eltype(TG))
end

function eliminatestaledates(lagged::AbstractVector{T}, date::AbstractVector{D},
  ::Type{T}, ::Type{D};
  laggeddate::AbstractVector{Union{D,Missing}}=error("laggeddate is required"),
  maxnotstale::Any=error("maxnotstale is required"),
  ) where {T<:Any, D<:Any}

  local missingindicies::Vector{Bool}

  if all((laggeddate .≤ date) .| (laggeddate .=== missing))
    missingindicies = ((ld, d)->
      ismissing(ld) || ismissing(d) || d - maxnotstale > ld).(laggeddate, date)
  #  println(DataFrame(laggeddate=laggeddate, date=date)[1:50,:])
  # with reversed dates, we are effectually computing a lead instead of a lag
  elseif all((laggeddate .≥ date) .| (laggeddate .=== missing))
    missingindicies = ((nd, d)->
      ismissing(nd) || ismissing(d) || nd - maxnotstale > d).(laggeddate, date)
    #@assert sum(missingindicies) ≠ length(lagged)
  else
    @assert false #something is wrong (df not sorted properly?)
  end

  lagged[missingindicies] .= missing

  #@assert !(all(missing .=== lagged))

  return lagged
end

function lagwithin2sorted(
  vals::AbstractVector,
  group::AbstractVector,
  date::AbstractVector{D};
  laggedgroup::Vector{Union{G, Missing}}=error("lagged group is required with lagwithin2sorted"),
  laggeddate::Vector{Union{D,Missing}}=error(
    "lagged group is required with lagwithin2sorted if maxnotstale is provided"),
  maxnotstale::Any = error(
    "maxnotstale is required if date vec is provided in lagwithin2sorted")) where {G<:Any, D<:Any}


  local lagged = lagwithin2sorted(vals, group, laggedgroup, eltype(vals), eltype(group))
  lagged = eliminatestaledates(lagged, date, eltype(lagged), eltype(date); laggeddate, maxnotstale)
  return lagged
end

#lags multiple columns within a dataframe
function lagwithin2sorted!(df::DataFrame, vals::Vector{<:DField}, group::DField;
  date::NDField = nothing,
  maxnotstale::Any = nothing,
  laggedvals::Vector{<:DField} = (s->Symbol(:L, s)).(vals))


  local laggedgroup::Vector{Union{eltype(df[!,group]), Missing}}
  local laggeddate::Vector{Union{eltype(df[!,date]), Missing}}


  #lag the groups, and dates if needed, once
  laggedgroup = Vector{Union{eltype(group), Missing}}([missing; df[!, group][1:(end-1)]])
  #allocate the taskes for multithreading
  tasks::Vector{Task} = Vector{Task}(undef, length(vals))

  #deploy the lags
  if maxnotstale === nothing
    tasks = @sync((s::DField -> Threads.@spawn lagwithin2sorted(
      df[!, s], df[!, group], laggedgroup, eltype(df[!, s]), eltype(df[!, group]),
      )).(vals))
  else
    laggeddate = [missing; df[!,date][1:(end-1)]]
    laggeddate[laggedgroup .!== df[!,group]] .= missing
    tasks = @sync((s::DField -> Threads.@spawn lagwithin2sorted(
      df[!, s], df[!, group], df[!,date]; laggedgroup, laggeddate, maxnotstale)).(vals))
  end


  for i ∈ 1:length(vals)
    df[!, laggedvals[i]] = fetch(tasks[i])
  end

  return nothing
end


function lagwithin2!(df::DataFrame, vals::Vector{<:DField}, group::DField;
  date::DField = error("date field is now required"),
  maxnotstale::Any = nothing,
  laggedvals::Vector{<:DField} = (s->Symbol(:L, s)).(vals))

  issorted(df, [group, date]) || error("df must be sorted by $group, $date")

  return lagwithin2sorted!(df, vals, group; date, maxnotstale, laggedvals)
end

function differencewithin2!(df::DataFrame,
  vals::Vector{<:DField},
  group::DField;
  date::DField = error("date field is now required"),
  differencedvals::Vector{<:DField} = (s::DField->Symbol(:D, s)).(vals),
  createlag::Bool=true,
  deletelag::Bool=true,
  maxnotstale::Any = nothing,
  laggedvals::Vector{<:DField} = (deletelag ?
    (s::DField->Symbol(:L, s, :_temp)).(vals) : (s->Symbol(:L, s)).(vals))
  )::Nothing

  (date ∈ vals) && error("date field cannot be differenced. workaround- create a copy")

  createlag && lagwithin2!(df, vals, group; date, laggedvals, maxnotstale)
  #println(df[1:50,:])
  for t ∈ 1:length(vals)
    df[!, differencedvals[t]] = df[!, vals[t]] .- df[!, laggedvals[t]]
  end
  deletelag && select!(df, Not(laggedvals))

  return nothing
end



function leadwithin2!(df::DataFrame, vals::Vector{<:DField}, group::DField;
  date::DField = error("date field is now required"),
  maxnotstale::Any = nothing,
  leadvals::Vector{<:DField} = (s->Symbol(:N, s)).(vals))

  issorted(df, [group, date]) || error("df must be sorted by $group, $date")


  #we reverse the target fields without copying
  inputfields = [date; vals; group;] |> unique!
  Nrows::Int = size(df,1)
  revdf = (f->f=>view(df[!,f], Nrows:-1:1)).(inputfields) |> (flds)->DataFrame(flds, copycols=false)
  #revdf = (f->f=>df[Nrows:-1:1,f]).(inputfields) |> DataFrame! (similar but copies columns)


  #now "lag" the reversed df (which is effectivelly a lead)
  lagwithin2sorted!(revdf, vals, group, date=date, maxnotstale=maxnotstale, laggedvals=leadvals)
  revrevdf = revdf[Nrows:-1:1,Not(inputfields)]
  for f ∈ propertynames(revrevdf)
    df[!, f] = revrevdf[!, f]
  end

  return nothing
end
