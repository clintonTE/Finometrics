


n2sStat(x::T where T<:Real) = num2Str(x, 3, Ints=true)

#calculates the test statistic and T value of a data vector
skewnessStat(data::Vector{Float64})::Float64 = skewness(data)/√(6. / length(data))
kurtosisStat(data::Vector{Float64})::Float64 = kurtosis(data)/√(24. / length(data))
JBStat(data::Vector{Float64})::Float64 =
  skewness(data)^2/(6. / length(data)) + kurtosis(data)^2/(24. / length(data))


#simple scheme for doing the ljungBox test
function ljungBox(ρ::Vector{Float64}, N::Int)::Float64
  m::Int = length(ρ)
  Q::Float64 = N*(N+2) * sum((ρ .* ρ) ./ (N .- collect(1:m)))
  return Q
end

#convenience function, defaults to log(T) for the lag
ljungBox(data::Vector{Float64}, lagRange::T where T <:AbstractRange = 1:(Int(round(log(length(data)))))) =
  ljungBox(autoρ(data, lagRange), length(data))


mutable struct UnivariateTest
  statistic::Function
  distribution::Distribution
  H₀::Float64
  testType::Symbol
end

#convenience constructor for univariate test
UnivariateTest(statistic::Function, distribution::Distribution;
  H₀::Float64=0.0, testType::Symbol=:upper)::UnivariateTest =
  UnivariateTest(statistic, distribution, H₀, testType)

function test(U::UnivariateTest, data::Vector{Float64}; verbose::Bool = true)::NTuple{2,Float64}
  θ::Float64 = U.statistic(data .- U.H₀)

  p::Float64 = cdf(U.distribution, θ)

  #adjust the p value for the specific test type
  if U.testType == :upper
    p = 1.0 - p
  elseif  U.testType == :lower
  elseif U.testType == :twoTailed
    p = p > 0.5 ? (1-p)*2 : 2*p
  else
    warn("Test type not found in univariate test")
  end


  return (θ, p)
end

#convenience function for test
test(U::UnivariateTest, data::T where T <: AbstractVector) = Vector{Float64}(dropmissing(data))


#simple function to get the autocorrelation
function autoρ(v::T, lag::Int)::Float64 where T <: AbstractVector
  vAvg::Float64 = mean(skipmissing(v))
  Δv::Vector{Float64} = v .- vAvg
  rng::S where S <:AbstractRange = (lag+1):length(Δv)
  vCov::Float64 = sum(((i::Int)->Δv[i]*Δv[i-lag]).(rng))
  return vCov/sum(Δv.^2.)
end

#convenience function to take in many lags
autoρ(v::T where T <: AbstractVector, lags::T where T<:AbstractRange)::Vector{Float64} = ((l::Int)->autoρ(v,l)).(lags)

#simple function to get the partial autocorrelation
function pAutoρ(v::T, lag::Int)::Float64 where T <: AbstractVector

  #fill the lag matrix
  Y::Vector{Float64} = v[(lag+1):end]
  X::Matrix{Float64} = ones(length(v)-lag,lag+1)
  for i::Int ∈ 1:lag
    X[:,i+1] .= v[(lag-i+1):(end-i)]
  end

  reg::FMLM = FMLM(X,Y)

  return reg.β[end]
end

#convenience function to take in many lags
pAutoρ(v::T where T <: AbstractVector, lags::T where T<:AbstractRange)::Vector{Float64} =
  ((l::Int)->pAutoρ(v,l)).(lags)

#Lag a series in a dataframe
function lagDF!(DF::DataFrame, target::Symbol, lags::Int = 1; name::Symbol = Symbol("$(target)L$lags"))::DataFrame
  DF[:,name] = similar(DF[:,target])
  DF[:,name] .= missing
  DF[(lags+1):end,name] = DF[1:(end-lags),target]

  return DF
end

function lagDF!(DF::DataFrame, targets::Vector{Symbol}, lags::Int = 1;
  names::Vector{Symbol}= ((s::Symbol)->Symbol("$(s)L$lags")).(targets))

  for i::Int ∈ 1:length(targets)
    lagDF!(DF, targets[i], lags, name = names[i])
  end

  return DF
end

mutable struct ARMA
  c::T where T <: Real
  ϕ::Vector{T}  where T <: Real
  ϑ::Vector{T} where T <: Real
  σ²::T where T <: Real
  dist::T where T<: Distribution
  p::Int
  q::Int
end

#constructor for ARMA
function ARMA(;c::Float64 = 0.0, ϕ::Vector{Float64}=ones(1), ϑ::Vector{Float64} = Vector{Float64}(),
    σ²::Float64 = 1.0, dist::UnivariateDistribution = Normal(0,σ²^0.5),
    p::Int = length(ϕ), q::Int = length(ϑ))::ARMA
  return ARMA(c,ϕ,ϑ,σ²,dist, p, q)
end

#this function allocates but does not paramaterize the model
function ARMA(p::Int, q::Int; c::Float64 = 0.0,
    σ²::Float64 = 1.0, dist::UnivariateDistribution = Normal(0,σ²^0.5))::ARMA

    return ARMA(c,zeros(p),zeros(q),σ²,dist, p, q)
end

#these are helper methods to change the type for things like autodifferentiation
function generalizeARMA!(Θ::ARMA)::Nothing
  Θ.c = Real(Θ.c)
  Θ.ϕ = Vector{Real}(Θ.ϕ)
  Θ.ϑ = Vector{Real}(Θ.ϑ)
  Θ.σ² = Real(Θ.σ²)

  return nothing
end

function narrowARMA!(Θ::ARMA)::Nothing
  Θ.c = Float64(Θ.c)
  Θ.ϕ = Vector{Float64}(Θ.ϕ)
  Θ.ϑ = Vector{Float64}(Θ.ϑ)
  Θ.σ² = Float64(Θ.σ²)

  return nothing
end

#fills the predicted Y values using the BJ methedology as described in Hamilton
function ARMAFillYHat!(Θ::ARMA, Y::Vector{Float64}, YHat::Vector{T} = similar(Y),
    ε::Vector{T} = similar(Y)
    #ε::Vector{T} = Vector{T}(length(Y)+max(Θ.q-Θ.p,0))
    )::NTuple{2,Vector{T}} where T<:Real

  ε .= 0.
  τ::Int = length(Y)

  #start with the observations to kick the model into gear
  YHat[1:Θ.p] .= Y[1:Θ.p]

  AR::T = 0.0
  MA::T = 0.0

  qIndexAdj::Int = max(Θ.q-Θ.p,0)

  #fill the predicted residual arrays
  for i::Int ∈ (Θ.p+1):τ
      for j::Int ∈ 1:Θ.p
          AR+=Y[i-j]*Θ.ϕ[j]
      end
      for j::Int ∈ 1:(min(Θ.q,i-1))
          MA -= ε[i-j]*Θ.ϑ[j]
      end
    #AR = Θ.p>0 ? sum(view(Y,(i-1):-1:(i-Θ.p)) .* Θ.ϕ):0.0
    #MA = Θ.q>0 ? sum( view(ε,(i-1+qIndexAdj):-1:(i-Θ.q+qIndexAdj)) .* Θ.ϑ):0.0
    #      println("got here")
    YHat[i] = AR + MA + Θ.c
    ε[i] = Y[i]-YHat[i]
    #ε[i+qIndexAdj] = Y[i] - YHat[i]
  end



  return (YHat, ε)
end



#this is from Hamilton
function ARMALogLike(Θ::ARMA, Y::Vector{Float64},YHat::Vector{T} = similar(Y),
    #ε::Vector{T} = zeros(T, length(Y)+max(Θ.q-Θ.p,0))
    ε::Vector{T} = similar(Y)
    )::T where T<:Real

  ARMAFillYHat!(Θ,Y,YHat,ε)

  #from foc of σ²
  ε²::T = sum((abs2).(ε))
  σ²::T = (ε²)/(length(Y)-Θ.p)


  #return -(length(Y)-Θ.p)/2*(log(2π)+log(Θ.σ²))-1/(2*Θ.σ²)*sum(abs2,ε[(Θ.p+1+max(Θ.q-Θ.p,0)):end])
  return -(length(Y)-Θ.p)/2*(log(2π)+log(σ²)) - ε²/(2*σ²)
end



function ARMAUpdate!(Θ::ARMA, θ::Vector{T}, suppressC::Bool)::ARMA where T<:Real
  if suppressC
    Θ.ϕ .= θ[1:(Θ.p)]
    Θ.ϑ .= θ[(Θ.p+1):end]
  else
    Θ.c = θ[1]
    Θ.ϕ .= θ[2:(Θ.p+1)]
    Θ.ϑ .= θ[(Θ.p+2):end]
  end

  return Θ
end


#holds estimate of an ARMA model
struct ARMAEstimate
  Θ::ARMA
  YHat::Vector{Float64}
  ε::Vector{Float64}
  cError::Float64
  ϕErrors::Vector{Float64}
  ϑErrors::Vector{Float64}
  σ²Error::Float64
end

#Algorithms: :LD_LBFGS , :LN_COBYLA, :LN_BOBYQA
#  :LD_SLSQPx , :LD_TNEWTON (fails), :LD_VAR2/:LD_VAR1 ,
#  :LD_MMA , :GN_ESCH, :GN_DIRECT , GN_DIRECT_L ,
# GN_ISRES, :LN_NELDERMEAD, :LN_SBPLX, :LD_LBFGS* (818) , :GD_STOGO
#GN_CRS2_LM
#estimates an ARMA model from the data
function ARMAEstimate(Θ::ARMA, Y::Vector{Float64};
    alg::Symbol = :LD_SLSQP, maxTime::Float64=MAX_TIME, fTolRel::Float64 = F_TOL_REL,
    fTolAbs::Float64 = F_TOL_ABS, bound::Float64 = OPT_BOUND, boundσ::Float64 = OPT_BOUND/10.,
    verbose::Bool = false, suppressC::Bool = false)::ARMAEstimate

  params::Vector{Float64} = suppressC ? [Θ.ϕ; Θ.ϑ] : [Θ.c; Θ.ϕ; Θ.ϑ]
  K::Int = length(params)

  #heuristics for the bound
  lowerBounds::Vector{Float64} = ones(K) .* -1. .* bound
  upperBounds::Vector{Float64} = ones(K) .* bound

  #safety check on parameters
  params .= max.(lowerBounds, (min.(upperBounds, params)))

  #preallocate by calling the fill YHat function, which will allocate the space
  YHat::Vector{Real}, ε::Vector{Real} = ARMAFillYHat!(Θ, Y)

  #use this for finite differencing
  Δz::Vector{Float64} = similar(params)

  #need to do this for autodifferentation
  generalizeARMA!(Θ)

  ##############objective functions
  function singleObj!(θ::Vector{T})::T where T<: Real
    Θ = ARMAUpdate!(Θ,θ, suppressC)
    return ARMALogLike(Θ,Y, YHat, ε)
  end

  function update∇!(θ::Vector{Float64}, ∇::Vector{Float64})::Nothing
    if length(∇) > 0
      ∇ .= ForwardDiff.gradient(singleObj!, θ) #∇Obj(θ)
    end
    return nothing
  end

  function obj!(θ::Vector{Float64}, ∇::Vector{Float64})::Float64
    objVal::Float64 = singleObj!(θ)
    update∇!(θ, ∇)

    return objVal
  end
  ###################optimization
  #println("testObj: $(obj!(params, similar(params)))")

  if verbose
    println("\nBeginning Optimization for ARMA($(Θ.p),$(Θ.q))")
    println("testObj: $(obj!(params, similar(params)))")
  end

  opt::NLopt.Opt = Opt(alg,K)

  max_objective!(opt, obj!)
  lower_bounds!(opt,lowerBounds)
  upper_bounds!(opt,upperBounds)
  #ftol_rel!(opt, fTolRel)
  #ftol_abs!(opt, fTolAbs)
  maxtime!(opt, maxTime)


  #run the actual optimziation
  (optVal,optθ,retCode) = optimize(opt, params)

  if verbose
    println("Objective: $optVal, retCode: $retCode")
    println("Parameters: $optθ")
  end
  #now do the standard errors
  #hessθ::Function = (θ::Vector{Float64})->ForwardDiff.hessian(singleObj!, θ)
  Σ::Matrix{Float64} = (-1 .* ForwardDiff.hessian(singleObj!, optθ)) \ eye(length(optθ)) #invert the hessian
  σ::Vector{Float64} = diag(Σ)
  if minimum(σ) < 0.
    println("WARNING: Negative Standard Variance Estimate, forced to 0")
  end
  σ .= ((x::Float64)->max(x,0.)^0.5).(σ)
  #println("σ2: $σ")

  if verbose

    Δh::Float64 = 10. ^ -9.
    ΔhMat::Matrix{Float64} = eye(length(optθ)) .* Δh
    function finite∇(θ::Vector{Float64}, ∇::Vector{Float64}=similar(θ), origObj::Float64 = singleObj!(θ))
      for i::Int=1:length(θ)
        ∇[i] = (origObj - singleObj!(θ .+ view(ΔhMat, :, i)))/(Δh)
      end

      return ∇
    end

    println("Analytical Gradient: $(ForwardDiff.gradient(singleObj!, optθ))")
    println("finite diff gradient: $(finite∇(optθ))")
  end

  narrowARMA!(Θ)
  ARMAUpdate!(Θ,optθ, suppressC)
  ARMAFillYHat!(Θ, Y, YHat, ε)
  Θ.σ² = sum((abs2).(Y .- YHat))/(length(Y)-Θ.p)
  Θ.dist = Normal(Θ.c, Θ.σ²^0.5)
  ΘEst::ARMAEstimate = suppressC ? ARMAEstimate(Θ, YHat, ε, 0.0, σ[1:Θ.p], σ[(Θ.p+1):end], 0.0) : (
    ARMAEstimate(Θ, YHat, ε, σ[1], σ[2:(Θ.p+1)], σ[(Θ.p+2):end], 0.0))

  return ΘEst
end


#quickly calculates the RMS of two vectors
RMS(v1::Vector{Float64}, v2::Vector{Float64})::Float64 =
  √(abs2(v1 .- v2)/length(v1))

#Based on StatsBase version, but handles missing in a different way
#Takes the prop percentile value, and replaces all lower values with that value
#Then takes the (1-prop) percentile value, and replaces all higher values with that value
#All calculations skip missing
function winsorize!(target::AbstractVector{T};
    prop::Float64 = 0.0, sorted::Bool = false)::Nothing where T <: Any

  local vcomplete::SubArray{T,1}

  if Missing <: T
    vcomplete =  view(target, (!ismissing).(target)) #don't trust how quantile handles missing
  else
    vcomplete=target
  end

  local minval::Float64 = quantile(vcomplete, min(prop,1-prop), sorted=sorted)
  local maxval::Float64 = quantile(vcomplete, max(prop,1-prop), sorted=sorted)

  target .= (f::MFloat64 -> max(minval, min(f,maxval))).(target)
  return nothing
end

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
