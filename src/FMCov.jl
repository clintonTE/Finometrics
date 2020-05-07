###################getHomosked!#######################
#Gets the covariance matrix under Homoskedastic assumptions σ^2[X'X]^-1

#main method to get the covariance matrix
#IN: The QR decomposition from the regression, a vector of residuals,
#optionally memory for the covariance matrix
#OUT: Writes and returns the covariance matrix
getHomoskedΣ!(args...; keyargs...) = error("use homoskedasticΣ instead")

function homoskedasticΣ!(xqr::FMQR{M}, ε::V, Σ::M = M(undef, lin.K, lin.K);
  dofcorrect::Float64 = 1.0)::M where  {M<:AbstractMatrix, V<:AbstractVector}

  BLAS.gemm!('N','T',(ε⋅ε)/xqr.N,xqr.Rinv,xqr.Rinv,0.0,Σ)

  Σ .*=

  return Σ #[X'X]^-1*σ2
end

#=helper method for the above function which extracts the QR decomposition
from a linear model
IN: A linear model, optionally memory for the covariance matrix
OUT: WRites and returns the covariance matrix=#
function homoskedasticΣ!(lin::FMLM{M},
  Σ::M = M(undef, lin.K, lin.K))::M where M<:AbstractMatrix

  return homoskedasticΣ!(lin.xqr, lin.ε, Σ, dofcorrect = lin.N / lin.dof)
end


#=Simple test function to verify the QR-related optimizations
Output should be identical to the standard methods
IN: Independent variable X matrix
OUT: Returns the covariance matrix=#
getHomoskedΣSlow(args...; keyargs...) = error("use homoskedasticΣslow instead")
function homoskedasticΣslow(X::M, ε::V; dofcorrect::Float64=1.0)::M where {M<:AbstractMatrix, V<:AbstractVector}
  T::Type = eltype(V)

  out = Matrix(X' * X)\I .* (ε⋅ε)/size(X,1)
  return out * dofcorrect
end

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
function homoskedasticΣslow(lin::FMLM{M, V}) where {M<:AbstractMatrix, V<:AbstractVector}
  return homoskedasticΣslow(lin.X, lin.Y.-lin.X * lin.β, dofcorrect = lin.N / lin.dof)::M
end


###################getWhiteERRORs ###OLS only
#checked against vcovHC in R (which uses the HC3 algorithm)
getWhiteΣ!(args...; keyargs...) = error("use whiteΣ! instead")
function whiteΣ!(xqr::FMQR{M}, ε::V,
  Σ::M = M(undef, xqr.K, xqr.K),
  dofcorrect::Float64=1.0)::M where {M<:AbstractMatrix, V<:AbstractVector}

  Λ::V = ε.^2.0

  QRtinv::M = BLAS.gemm('N','T',xqr.Q,xqr.Rinv) #Q(R')^-1
  RinvQΛ::M = QRtinv'
  RinvQΛ *= Diagonal(Λ)

  BLAS.gemm!('N','N',1.0,RinvQΛ, QRtinv, 0.0, Σ)

  Σ .*= dofcorrect

  #final multiplicaiton and assignment
  return Σ #R^-1Q'ΛQ(R')^-1
end

function whiteΣ!(lin::FMLM{M},
  Σ::M = M(undef, lin.K, lin.K))::M where M<:AbstractMatrix

  return whiteΣ!(lin.xqr, lin.ε, Σ, lin.N/lin.dof)
end

#for testing purposes only
whiteΣslow(args...; keyargs...) = error("use whiteΣslow instead")
function whiteΣslow(X::M, ε::V, dofcorrect::Float64=1.0)  where {M<:AbstractMatrix, V<:AbstractVector}
  XXinv::M = (X' * X)\I

  return XXinv * X' * diagm(ε) .^ 2.0 * X * XXinv .*dofcorrect
end

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
whiteΣslow(lin::FMLM) = whiteΣslow(lin.X, lin.Y.-lin.X * lin.β, lin.N/lin.dof)

###########################getClustered SEs
#method for getting clustered errors
#from http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf
#also provides methedology for multi-way clsutering
#forumala: [X'X]^-1*B*[X'X]^-1 where B= sum over G (Xg'*εg*εg'*Xg) and g is indexed for G clusters
#checked against felm's [modelvar].cse with exactDOF=TRUE in R
getClustered!(args...; keyargs...) = error("use clusteredΣ! instead")

function clusteredΣ!(X::M, xqr::FMQR{M}, ε::V, clusters::C,
    Σ::M = M(undef, xqr.K, xqr.K),
    dof::Int = xqr.N-xqr.K)::M where {M<:AbstractMatrix, V<:AbstractVector, C<:AbstractVector}

  #convenience vectors and values
  inddf::DataFrame = DataFrame(clusters=clusters, idx = 1:xqr.N)
  sindfs::GroupedDataFrame = groupby(inddf, :clusters)
  Nclusters::Int = length(sindfs)

  #ctable::Dict = Dict(cvars[i] => i for i::Int ∈ 1:Nclusters)
  #clustercode::Vector{Int} = ((x)->ctable[x]).(clusters)

  #iterate through the groups
  B::M = zeros(xqr.K, xqr.K)
  @fastmath for sindf::SubDataFrame ∈ sindfs
    Xg::SubArray = view(X, sindf.idx, :)
    ug::SubArray = view(ε, sindf.idx)
    B .+= Xg'*ug*ug'*Xg
  end

  # Calc [X'X]^-1*B*[X'X]^-1
  RRinv::M =  BLAS.gemm('N','T',xqr.Rinv,xqr.Rinv)
  RRinvB::M = BLAS.gemm('N','N',RRinv,B)

  dofcorrect::Float64 = Nclusters/(Nclusters-1.)*(xqr.N-1.)/dof #small sample correction

  if M<:Matrix{Float64}
    BLAS.gemm!('N','N',dofcorrect,RRinvB,RRinv,0.0,Σ)
  else
    Σ = Matrix(BLAS.gemm('N','N',dofcorrect,RRinvB,RRinv))
  end

  return Σ
end

#handles 2-way clustering
#checked against felm's [modelvar].cse with exactDOF=TRUE in R
function clusteredΣ!(lin::FMLM{M, V}, Σ₁::M = M(undef, lin.K, lin.K);
    clusters::Vector{<:Vector}#=Union{FMClusters, Vector{<:FMData}}=# = lin.clusters, #allows for an override on the standard errors
    testequivelance::Bool = false)::M where {M<:AbstractMatrix, V<:AbstractVector}


  length(clusters) ∈ [1,2] || error("Illegal number of clusters = $(length(clusters))")
  clusteredΣ!(lin.X, lin.xqr, lin.ε, clusters[1], Σ₁, lin.dof)
  if length(clusters) == 2
    Σ₂::M = similar(Σ₁)
    Σ₃::M = similar(Σ₁)

    #create the intersection of the clusters

    cluster⋂::Vector{Symbol} = ((s1, s2)->
      Symbol(string(s1),"_∩_",string(s2))).(clusters[1],clusters[2])

    #get the remaining clusters
    clusteredΣ!(lin.X, lin.xqr, lin.ε, clusters[2], Σ₂, lin.dof)

    #may verify the equivelance of white Σ and the cluster intersection Σ in the case where
    #cluster 1 within cluster 2 accounts for all variation (e.g. combination of clusters is unique)
    if testequivelance
      clusteredΣ!(lin.X, lin.xqr, lin.ε, cluster⋂, Σ₃, lin.dof)
      if length(unique(cluster⋂)) ≠ length(cluster⋂)
        @warn "Note: cluster intersection not unique. White results will not reconcile."
      end
      println("Clustered version:")
      display(Σ₃)

      whiteΣ!(lin.xqr, lin.ε, Σ₃, lin.N/lin.dof)
      println("White version:")
      display(Σ₃)
    else #standard case
      if length(unique(cluster⋂)) == length(cluster⋂) #special case allowing for performance optimization
        whiteΣ!(lin.xqr, lin.ε, Σ₃, lin.N/lin.dof)
      else
        clusteredΣ!(lin.X, lin.xqr, lin.ε, cluster⋂, Σ₃, lin.dof)
      end
    end

    Σ₁ .= Σ₁ .+ Σ₂ .- Σ₃

  end

  return Σ₁
end



###################getMWErrors!##############################
#Main method to get the modified white SEs Does #[X'X]^-1X'ΛX[X'X]^-1


#=Main method to calculate the modified white SEs
IN: A FMQR decomposition object, a residuals vector,
optional memory for a covariance matrix
OUT: Writes and returns the covariance matrix=#
getModWhiteΣ!(args...; keyargs...) = error("Use modifiedwhiteΣ! instead")
function modifiedwhiteΣ!(xqr::FMQR{M}, ε::V,
  Σ::M = M(undef, xqr.K, xqr.K), dofcorrect::Float64 = 1.0
  )::M where {M<:AbstractMatrix, V<:AbstractVector}

  Λ::V = similar(ε)

  project!(xqr,Λ) #get the projection diagonal
  Λ .= (ε./(1.0.-Λ)).^2.0

  QRtinv::M = BLAS.gemm('N','T',xqr.Q,xqr.Rinv) #Q(R')^-1
  RinvQΛ::M = QRtinv'

  #loop to scale the matrix by Λ (modified part of modified white)
  #=@fastmath for j ∈ 1:xqr.N, i∈1:xqr.K
    RinvQΛ[i,j] *= Λ[j]
  end=#
  RinvQΛ *= Diagonal(Λ) #hopefully this is fast
  T::Type = eltype(M)

  BLAS.gemm!('N','N',dofcorrect,RinvQΛ, QRtinv,0.0,Σ)
  #final multiplicaiton and assignment
  return  Σ
end

#=helper method which extracts the required components from a linear model
IN: A linear model, optionally memory for the covariance matrix
OUT: Writes and allocates the covariance matrix=#
function modifiedwhiteΣ!(lin::FMLM{M},
  Σ::M = M(undef, lin.K, lin.K))::M where M <: AbstractMatrix

  return modifiedwhiteΣ!(lin.xqr, lin.ε, Σ, lin.N/lin.dof)
end


#=Simple test function to verify the QR-related optimizations
Output should be identical to the standard methods
IN: Independent variable X matrix
OUT: Returns the covariance matrix=#
getModWhiteΣSlow(args...; keyargs...) = error("Use modifiedwhiteΣslowinstead")
function modifiedwhiteΣslow(X::M, ε::V, dofcorrect::Float64) where {M<:AbstractMatrix, V<:AbstractVector}

  #do this since we only care about clarity (and not performance)
  (M<:CuArray) && CuArrays.allowscalar(true)

  XXinv::M = Matrix(X' * X)\I

  P::M = X * (XXinv) * X'
  Σ::M = XXinv * X' * diagm(ε  ./ (1.0 .- diag(P))) .^ 2.0 * X * XXinv

  (M<:CuArray) && CuArrays.allowscalar(false)


  return Σ .* dofcorrect

end

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
modifiedwhiteΣslow(lin::FMLM, dofcorrect::Float64 = lin.N/lin.dof) = modifiedwhiteΣslow(
  lin.X, lin.Y.-lin.X * lin.β, dofcorrect)

  #########################getNeweyWest ##########OLS Only
  #NOTE: To match Sandwich NeweyWest in R, set prewhite = FALSE, adjust = TRUE
  getNeweyWest!(args...; keyargs...) = error("use neweywestΣ! instead")

  function neweywestΣ!(X::M, xqr::FMQR{M}, ε::V, lag::Int,
    Σ::M = M(undef, xqr.K, xqr.K),
    dofcorrect::Float64 = 1.0)::M where {M<:AbstractMatrix, V<:AbstractVector}

    #pre-allocate for the spectral matrix

    #NOTE: we need to switch away from gpu arrays due to the lagging
    #that is, views don't seem to work properly on gpu matrices
    T::Type = eltype(M)
    Rv::Matrix{T} = Matrix{T}(undef, xqr.K, xqr.K) #pre-allocate working matrix
    RRinv::Matrix{T} = BLAS.gemm('N', 'T', xqr.Rinv, xqr.Rinv) #this is equivelent to [X'X]^-1

    #need to multiply through by the error
    Xe::Matrix{T} = X .* ε
    ST::Matrix{T} = BLAS.gemm('T','N',T(1.0/xqr.N), Xe, Xe)
    for v::Int ∈ 1:lag
      #overwrites Rv with (1/N)R'R
      BLAS.gemm!('T', 'N', T(1.0/xqr.N), view(Xe, (v+1):(xqr.N), :),view(Xe, 1:(xqr.N-v), :), T(0.0), Rv)
      #Rv .= view(Xe, (v+1):(xqr.N), :)' * view(Xe, 1:(xqr.N-v), :) .* (1.0/xqr.N)
      ST += (lag + 1 - v)/(lag+1.) .* (Rv .+ Rv')
    end

    #this is [X'X]^-1S=[R'R]^-1S
    RRinvS::Matrix{T} = BLAS.gemm('N', 'N', RRinv, ST)

    #finally we have T[X'X]^-1S[X'X]^-1
    if M<:Matrix{Float64}
      BLAS.gemm!('N','N',Float64(xqr.N), RRinvS, RRinv, 0.0, Σ)
    else
      Σ = BLAS.gemm('N','N',T(xqr.N), RRinvS, RRinv)
    end


    Σ .*= dofcorrect
    #println("$(diag(Σ .* dofcorrect))")sasa
    return Σ
  end


  function neweywestΣ!(lin::FMLM{M,V}, lag::Int,
    Σ::M = M(undef, lin.K, lin.K))::M where {M<: AbstractMatrix, V<:AbstractVector}

    return neweywestΣ!(lin.X, lin.xqr, lin.ε, lag, Σ, lin.N/lin.dof)
  end

  #helper function in case the input requires a single argument function
  getNeweyWestFunc(args...; keyargs...) = error("use neweywestΣfunc instead")
  neweywestΣfunc(lag::Int) = (lin::FMLM)-> neweywestΣ!(lin, lag)

  #testing for NeweyWest
  #NOTE: assumes dofcorrect value
  getNeweyWestSlow(args...; keyargs...) = error("use neweywestΣslow instead")
  function neweywestΣslow(X::M, ε::V, lag::Int, dofcorrect::Float64
      )  where {M<:AbstractMatrix, V<:AbstractVector}

    N::Int, K::Int = size(X)
    (M<:CuArray) && CuArrays.allowscalar(true) #don't care about performance here

    #first form the spectral matrix

    #initial value for spectral matrix
    xxt::Matrix = zeros(K,K)
    for t::Int ∈ 1:N
      xxt .+= (X[t, :] * X[t, :]') .* ε[t]^2
    end
    ST::Matrix =  xxt ./ N

    #make the rest of the spectral matrix
    for v::Int ∈ 1:lag
      xxt .= zeros(K,K)
      for t::Int ∈ (1+v):(N)
        xxt .+= (X[t, :] * X[t-v, :]') .* ε[t]*ε[t-v]
      end
      xxt ./= N
      ST .+= (lag + 1. - v)/(lag + 1.) .* (xxt .+ xxt')
    end

    XXinv::Matrix = Matrix(X' * X) \ I

    ret::M = N*XXinv*ST* XXinv
    (M<:CuArray) && CuArrays.allowscalar(false) #don't care about performance here
    return ret * dofcorrect
  end



  neweywestΣslow(lin::FMLM, lag::Int) = neweywestΣslow(lin.X, lin.Y .- lin.X * lin.β, lag, lin.N/lin.dof)

  getNeweyWestFuncSlow(args...; keyargs...) = error("use neweywestΣslow instead")
  neweywestΣfuncslow(lag::Int)::Function = (lin::FMLM)-> neweywestΣslow(lin, lag)

  ######################Panel Newy-West#######################

function neweywestpanelΣ!(X::M, xqr::FMQR{M}, ε::V, clusters::FMClusters,
    lag::Int, Σ::M = M(undef, xqr.K, xqr.K))::M where {
    M<:AbstractMatrix, V<:AbstractVector}

  if length(clusters)==0
    error("No clusters detected in neweywestpanelΣ! - use neweywestΣ! instead")
  elseif length(clusters)>1
    error("Multi-way clusters not supported in neweywestpanelΣ!")
  end

  N::Int = length(ε)
  E::Type = eltype(M)

  K::Int = size(X,2)
  Sₜ::M = zeros(K,K) #holds teh central matrix

  temp::M = M(undef, K, K) #pre-allocate working matrix
  RRinv::M = BLAS.gemm('N', 'T', xqr.Rinv, xqr.Rinv) #this is equivelent to [X'X]^-1

  inddf::DataFrame = DataFrame(clusters=clusters[1], idx = 1:xqr.N)
  sindfs::GroupedDataFrame = groupby(inddf, :clusters)
  Nclusters::Int = length(sindfs)

  #iterate over all stocks
  for sindf::SubDataFrame ∈ sindfs
    Xₙ::SubArray = view(X, sindf.idx, :)
    εₙ::SubArray = view(ε, sindf.idx)

    #need to multiply through by the error
    T::Int = length(εₙ)
    Xₑ::M = Xₙ .* εₙ

    Sₜ .+= BLAS.gemm('T','N',  E(1.0/N), Xₑ, Xₑ)

    for v::Int ∈ 1:lag
      #overwrites temp with (1/N)R'R
      BLAS.gemm!('T', 'N',   E(1.0/N), view(Xₑ, (v+1):T, :),view(Xₑ, 1:(T-v), :),   E(0.0), temp)
      Sₜ .+= (lag + 1 - v)/(lag+1.) .* (temp .+ temp')
    end

  end
  #Sₜ = Matrix{Float64}(I,K,K)
  #this is [X'X]^-1S=[R'R]^-1S
  RRinvS::M = BLAS.gemm('N', 'N', RRinv, Sₜ)

  #display(RRinv)

  #finally we have T[X'X]^-1S[X'X]^-1
  BLAS.gemm!('N','N', E(N), RRinvS, RRinv, E(0.0), Σ)
  #println("$(diag(Σ .* dofCorrect))")sasa
  return Σ
end

function neweywestpanelΣ!(lin::FMLM{M,V}, lag::Int,
  Σ::M = M(undef, lin.K, lin.K))::M where {M<: AbstractMatrix, V<:AbstractVector}

  return neweywestpanelΣ!(lin.X, lin.xqr, lin.ε, lin.clusters, lag, Σ)
end


neweywestpanelΣfunc(lag::Int) = (lin::FMLM)-> neweywestpanelΣ!(lin, lag)


function neweywestpanelΣslow!(X::M, xqr::FMQR{M}, ε::V, clusters::FMClusters,
    lag::Int, Σ::M = M(undef, xqr.K, xqr.K))::M where {
    M<:AbstractMatrix, V<:AbstractVector}

  N::Int = size(X,1)
  K::Int = size(X,2)
  T::Type = eltype(M)

  if length(clusters)==0
    error("No clusters detected in neweywestpanelΣslow! - use neweywestΣslow! instead")
  elseif length(clusters)>1
    error("Multi-way clusters not supported in neweywestpanelΣslow!")
  end

  #pre-allocate
  Sₜ::M = zeros(K, K)

  #this method requires a sorted array
  #will do this via dataframes
  df::DataFrame = DataFrame(deepcopy(X))
  xnames = names(df)
  df.clusters = deepcopy(clusters[1])
  df.ε = deepcopy(ε)

  sort!(df, :clusters)
  X = M(df[!, xnames])
  ε = V(df.ε)
  clusters = df.clusters

  #helper function using the Lars method
  for t ∈ 1:N
    Sₜ .+= X[t,:] * X[t,:]' * ε[t]^2
  end

  for v::Int ∈ 1:lag
    for t ∈ (v+1):N
      (clusters[t] == clusters[t-v]) && (
        Sₜ .+= (lag + 1 - v)/(lag+1.) .* (X[t,:] * X[t-v,:]' .+ X[t-v,:] * X[t,:]') * ε[t] * ε[t-v])
    end
  end
  Sₜ ./= N
  #Sₜ = Matrix{Float64}(I,K,K)

  #the final step is to multiply out the var-covar sandwhich
  RRinv::M = BLAS.gemm('N', 'T', xqr.Rinv, xqr.Rinv) #this is equivelent to [X'X]^-1

  #this is [X'X]^-1S=[R'R]^-1S
  RRinvS::M= BLAS.gemm('N', 'N', RRinv, Sₜ)

  #finally we have T[X'X]^-1S[X'X]^-1
  BLAS.gemm!('N','N',T(N), RRinvS, RRinv, T(0.0), Σ)

  return Σ
end

function neweywestpanelΣslow!(lin::FMLM{M,V}, lag::Int,
  Σ::M = M(undef, lin.K, lin.K))::M where {M<: AbstractMatrix, V<:AbstractVector}

  return neweywestpanelΣslow!(lin.X, lin.xqr, lin.ε, lin.clusters, lag, Σ)
end

neweywestΣslowfunc(lag::Int) = (lin::FMLM)-> neweywestpanelΣslow!(lin, lag)
