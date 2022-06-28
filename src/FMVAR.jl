


#holds prediction information
mutable struct VARPrediction
  dates::Vector{T} where T<:Union{Date, Int, Float64}
  X::Matrix{Float64}
  Y::Matrix{Float64}
  YHat::Matrix{Float64}
end

#holds fundamental information on the VAR
mutable struct VAR
  ϕ₀::Vector{Float64}
  ϕ₁::Matrix{Float64}
  compactϕ₀::Vector{Float64}
  compactϕ₁::Matrix{Float64}
  μ::Vector{Float64}
  L::Matrix{Float64}
  structures::Vector{T}  where T <: Structure
end

#holds the VAR and any estimation information
mutable struct VAREstimate
  var::VAR
  coeffErrors::Vector{T} where T <: CoefficientErrors
  N::Vector{T} where T <: RegressionN
  SSE::Vector{T} where T <: RegressionSS
  SSR::Vector{T} where T <: RegressionSS
  SST::Vector{T} where T <: RegressionSS
end


function VAR(ϕ₀::Vector{Float64}, ϕ₁::Matrix{Float64}, L::Matrix{Float64},
    structures::Vector{T}  where T <: Structure)::VAR

  #calculate the mean
  K::Int = length(ϕ₀)
  μ::Vector{Float64} = ((Matrix{Float64}(I,K,K) .- ϕ₁)\Matrix{Float64}(I,K,K))*ϕ₀

  compactRows::BitArray = structures .== nothing

  compactϕ₀::Vector{Float64} = ϕ₀[compactRows]
  compactϕ₁::Matrix{Float64} = ϕ₁[compactRows,compactRows]

  return VAR(ϕ₀, ϕ₁, compactϕ₀, compactϕ₁, μ,  L, structures)
end


#not sure how well this works with restrictions, if at all
function propagateImpulse(var::VAR, Z::Vector{Float64}, lags::Int=20)::Matrix{Float64}
  compactShocks::Vector{Float64} = var.L * Z
  K::Int = length(var.ϕ₀)
  propagate::Matrix{Float64} = zeros(K, lags + 1)

  #expand the shock vector
  compactRows::BitArray = var.structures .== nothing
  propagate[compactRows,1] .= compactShocks

  #this accounts for the demeaning

  for i::Int ∈ 1:lags
    propagate[:, i+1] = var.ϕ₁*propagate[:,i]
  end

  return propagate
end

#make period or multi-period ahead perdictions
function VARPrediction(var::VAR, dates::Vector{D}, X::Matrix{Float64}, Y::Matrix{Float64};
    τ::Int = 1)::VARPrediction where D <: Union{Date, Int, Float64}

  K::Int = length(var.ϕ₀)
  T::Int = size(Y,2)

  #housekeeping with lengths
  Y=Y[:,(τ+1):end]
  dates = dates[(τ+1):end]
  YHat::Matrix{Float64} = Matrix{Float64}(K, T-τ)
  ε::Matrix{Float64} = similar(YHat)
  predict::Vector{Float64} = Vector{Float64}(K)

  for i::Int ∈ 1:(T-τ)
    predict .= var.ϕ₀ .+ var.ϕ₁ * X[:,i] #first period ahead is based on data

    for j::Int ∈ 2:τ #subsequent periods based on previous prediction
      predict .= var.ϕ₀ .+ var.ϕ₁ * predict
    end
    YHat[:,i] .= predict
    ε[:,i] = Y[:, i] .- predict
  end

  return VARPrediction(dates, X, Y, YHat)
end


function varΣ(var::VAR)::Matrix{Float64}
  K::Int = length(var.compactϕ₀)

  #maybe a more efficient way to do this?
  bigMat::Matrix{Float64} = (eye(K^2) .- kron(var.compactϕ₁, var.compactϕ₁))\eye(K^2)
  Γ₀::Matrix{Float64} = reshape(bigMat*vec(var.L*var.L'),K,:)

  return Γ₀
end


#this is the low-level VAR estimator
function VAREstimate(DF::DataFrame, LHS::Vector{Symbol}, RHS::Vector{Symbol};
      errors::Function = getHomoskedΣ!, K::Int = length(LHS),
      structures::Vector{Structure}  = Vector{Structure}(fill(nothing, K)),
      ϕ::Matrix{Float64}=Matrix{Float64}(undef, K,K+1), #these are available for pre-allocaiton purposes
      σ::Vector{Float64} = Vector{Float64}(undef, K))

  #build up a couple of indices and helper variables
  K::Int = length(LHS)
  RHSIndex::Dict = Dict(RHS[i]=>i for i::Int ∈ 1:K)
  LHSIndex::Dict = Dict(LHS[i]=>i for i::Int ∈ 1:K)

  coeffErrors::Vector{CoefficientErrors} = Vector{Structure}(fill(nothing, K))
  N::Vector{RegressionN} = Vector{RegressionN}(fill(nothing, K))
  ε::Vector{Vector{Float64}} = Vector{Vector{Float64}}()

  #store the ANOVA information
  SSE::Vector{RegressionSS} = Vector{RegressionSS}(fill(nothing, K))
  SSR::Vector{RegressionSS} = similar(SSE)
  SST::Vector{RegressionSS} = similar(SSE)

  XExpr::CTExpr = Meta.parse(vec2String(RHS, " + "))

  #check for structure. If its not a structure (either specified or implied as a lag),
  #run a regression and estimate the coefficients
  for k::Int ∈ 1:K
    if structures[k] ≠ nothing
      ϕ[k,:] .= structures[k]
      σ[k] = 0.0
      coeffErrors[k], N[k], ε[k], SSE[k], SSR[k], SST[k]= fill(nothing, 6)
    elseif haskey(RHSIndex,LHS[k]) #by defintion, the lag must be restricted
      ϕ[k,:] .= 0.0
      ϕ[k,RHSIndex[LHS[k]]+1] = 1.0 #lookup the location and restict the outcome to an identity
      structures[k] = Vector{Float64}(ϕ[k,:])
      σ[k] = 0.0
      coeffErrors[k], N[k], ε[k], SSE[k], SSR[k], SST[k] = fill(nothing, 6)
    else
      structures[k] = nothing #probably redundant
      reg::CTLM = CTLM(DF,  XExpr, LHS[k], XNames=[:Intercept; RHS], YName=LHS[k])
      ϕ[k,:] .= reg.β
      push!(ε, reg.ε)
      σ[k] = std(ε[end])
      coeffErrors[k] = diag(getHomoskedΣ!(reg)).^0.5
      N[k] = reg.N
      SSE[k] = sum(reg.ε.^2)
      SSR[k] = sum((reg.Y .- reg.ε .- mean(reg.Y)).^2)
      SST[k] = sum((reg.Y .- mean(reg.Y)).^2)
    end
  end

  L::Matrix{Float64} = (cholesky(cov(hcat(ε...)))).L #(splat and calculate the cholesky root)

  var::VAR = VAR(Vector{Float64}(ϕ[:,1]), ϕ[:,2:end], L, structures)

  #note this only works for the linearly independent components



  return   VAREstimate(var, coeffErrors, N, SSE, SSR, SST)
end
