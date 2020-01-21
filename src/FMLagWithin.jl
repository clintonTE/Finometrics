
#lags a variable in the time series
#assumes the order field is intended to be in ascending order
function lagwithin!(df::DataFrame,
  targets::Vector{Symbol},
  groups::Vector{Symbol},
  period::Symbol;
  lags::Int=1,
  laggednames::Vector{Symbol} = (lags ≠ 1 ?
    (s::Symbol->Symbol(:L, lags, s)).(targets) : (s::Symbol->Symbol(:L, s)).(targets)),
  sorted::Bool=false)::Nothing

  T::Type = eltype(df[!,period]) #the type of the period column

  if !sorted
    sort!(df, [groups; period])
  end

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
  targets::Union{Symbol, Vector{Symbol}},
  groups::Union{Symbol, Vector{Symbol}},
  period::Symbol;
  lags::Int=1,
  laggedname::NSymbol = nothing,
  laggednames::Vector{Symbol} = (lags ≠ 1 ?
    (s::Symbol->Symbol(:L, lags, s)).([targets;]) : (s::Symbol->Symbol(:L, s)).([targets;])),
  sorted::Bool=false)::Nothing = (lagwithin!(df, [targets;], [groups;], period,
      lags=lags, laggednames=[something(laggedname, laggednames);], sorted=sorted))


#creates a differenced column
function differencewithin!(df::DataFrame,
  targets::Vector{Symbol},
  groups::Vector{Symbol},
  period::Symbol;
  differencednames::Vector{Symbol} = (s::Symbol->Symbol(:D, s)).(targets),
  sorted::Bool=false,
  createlag::Bool=true,
  deletelag::Bool=true,
  laggednames::Vector{Symbol} = (deletelag ?
    (s::Symbol->Symbol(:L, s, :_temp)).(targets) : (s->Symbol(:L, s)).(targets))
  )::Nothing

  createlag && lagwithin!(df, targets, groups, period,
    laggednames=laggednames, sorted=sorted)
  for t ∈ 1:length(targets)
    df[!, differencednames[t]] = df[!, targets[t]] .- df[!, laggednames[t]]
    deletelag && select!(df, Not(laggednames[t]))
  end

  return nothing
end

#helper method to handle the case of a single target and/or group
differencewithin!(df::DataFrame,
  targets::Union{Symbol, Vector{Symbol}},
  groups::Union{Symbol, Vector{Symbol}},
  period::Symbol;
  differencedname::NSymbol = nothing,
  differencednames::Vector{Symbol} = (s::Symbol->Symbol(:D, s)).([targets;]),
  sorted::Bool=false,
  createlag::Bool=true,
  deletelag::Bool=true,
  laggedname::NSymbol = nothing,
  laggednames::Vector{Symbol} = (deletelag ?
    (s::Symbol->Symbol(:L, s, :_temp)).([targets;]) : (s->Symbol(:L, s)).([targets;]))
  )::Nothing = differencewithin!(df, [targets;], [groups;], period,
    sorted=sorted,  createlag=createlag, deletelag=deletelag,
    differencednames=[something(differencedname, differencednames);],
    laggednames=[something(laggedname, laggednames);])


##################Updated versions

###if these work, consider moving to finometrics
#lag within a group
function lagwithin2sorted2(vals::AbstractVector{T},
  group::AbstractVector{G}, laggedgroup::Vector{<:Union{G, Missing}},
  ::Type{T}, ::Type{G}) where {T<:Any, G<:Any}

  lagged::Vector{Union{T, Missing}} = [missing; vals[1:(end-1)]]
  missingindicies::Vector{Bool} = ((lg,g)->ismissing(lg) || ismissing(g) || lg≠g).(laggedgroup,group)
  lagged[missingindicies] .= missing

  return lagged
end

function eliminatestaledates(lagged::AbstractVector{T}, date::AbstractVector{D},
  laggeddate::AbstractVector{Union{D,Missing}}, maxnotstale::Any,
  ::Type{T}, ::Type{D}) where {T<:Any, D<:Any}

  missingindicies::Vector{Bool} = ((ld, d)->
    ismissing(ld) || ismissing(d) || d - ld > maxnotstale).(laggeddate, date)
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


  local lagged = lagwithin2sorted2(vals, group, laggedgroup, eltype(vals), eltype(group))

  if !isnothing(maxnotstale)
    if isnothing(laggeddate)
      laggeddate = [missing; [date[1:(end-1)]]]
    end
    lagged = eliminatestaledates(lagged, date, laggeddate, maxnotstale, eltype(lagged), eltype(date))
  end

  return lagged
end

#lags multiple columns within a dataframe
function lagwithin2sorted!(df::DataFrame, vals::Vector{Symbol}, group::Symbol;
  date::Union{Symbol,Nothing} = nothing,
  maxnotstale::Any = nothing,
  laggedvals::Vector{Symbol} = (s->Symbol(:L, s)).(vals))

  local laggedgroup::Vector{Union{eltype(df[!,group]), Missing}}
  local laggeddate::Union{Vector{Union{eltype(df[!,date]), Missing}}, Nothing}


  #lag the groups, and dates if needed, once
  laggedgroup = Vector{Union{eltype(group), Missing}}([missing; df[!, group][1:(end-1)]])
  if !isnothing(date)
    laggeddate = [missing; df[!,date][1:(end-1)]]
  else
    laggeddate = nothing
  end
  #println(typeof(laggeddate))
  #allocate the taskes for multithreading
  tasks::Vector{Task} = Vector{Task}(undef, length(vals))

  #println(typeof(laggeddate))
  #deploy the lags
  tasks = (s::Symbol -> Threads.@spawn lagwithin2sorted(
    df[!, s], df[!, group], date = isnothing(date) ? nothing : df[!,date],
    laggeddate = laggeddate, maxnotstale = maxnotstale, laggedgroup = laggedgroup)).(vals)
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


function lagwithin2!(df::DataFrame, vals::Vector{Symbol}, group::Symbol;
  date::Union{Symbol,Nothing} = nothing,
  maxnotstale::Any = nothing,
  laggedvals::Vector{Symbol} = (s->Symbol(:L, s)).(vals),
  sorted::Bool = false)

  (!sorted) && sort!(df, [group, date])

  return lagwithin2sorted!(df, vals, group, date=date, maxnotstale=maxnotstale, laggedvals=laggedvals)
end



function testlagwithin2(N=100_000)

  #make the test data
  df::DataFrame = DataFrame(
    val1 = Vector{MFloat64}(undef, N),
    val2 = Vector{MFloat64}(undef, N),
    group = rand(collect(1:1000),N),
    date = rand(collect(Date(1985,11,11):Day(1):Date(2011,11,11))))

  for i ∈ 1:N, s ∈ [:val1, :val2]
    if rand() < 0.95
      df[i, s] = missing
    else
      df[i, s] = rand()
    end
  end

  maxnotstale = Day(1000)

  #execute! (will also sort)
  lagwithin2!(df, [:val1, :val2], :group, date=:date, sorted=false, maxnotstale = maxnotstale)

  df.TLval1 = similar(df.Lval1)
  df.TLval2 = similar(df.Lval2)
  #make the test lags
  for sdf ∈ groupby(df, :group)
    for j ∈ 2:size(sdf,1)
      if sdf[j-1,:date] ≥ sdf[j,:date] - maxnotstale
        sdf[j,:TLval1] = sdf[j-1,:val1]
        sdf[j,:TLval2] = sdf[j-1,:val2]
      end
    end
  end

  inequalities::Int = 0
  for r ∈ eachrow(df)
    if ismissing(r.Lval1)
      inequalities += (!ismissing(r.TLval1))
    else
      inequalities += (!(r.Lval1==r.TLval1))
    end

    if ismissing(r.Lval2)
      inequalities += (!ismissing(r.TLval2))
    else
      inequalities += (!(r.Lval2==r.TLval2))
    end
  end

  @assert inequalities == 0

  @info "passed lagwithin2 test"
end
