
##################Utility Functions ###############################

#helper function  to get the model matrix from a dataframe given a formula
#IN: A dataframe and formula object
#OUT: the model matrix
function generateX(::Type{M}, df::T, f::FormulaTerm)::M where
  {M<:AbstractMatrix, T <: AbstractDataFrame}

  #*************WARNING WARNING WARNING WARNING WARNING
  #replace the code below w/ the two lines above
  #once the missing-value coding of numerical variables is fixed see
  #https://github.com/JuliaStats/StatsModels.jl/issues/145

  #f = apply_schema(f, schema(f,df), StatisticalModel)
  #m = modelcols(f.rhs, df)
  (M<:CuArray) && CuArrays.allowscalar(true)
  m::M = ModelMatrix(ModelFrame(f, df)).m
  (M<:CuArray) && CuArrays.allowscalar(false)

  return m
end

#Same as above but allows for an expression
function generateX(::Type{M}, df::T, rhs::V)::M where
  {M<:AbstractMatrix, T <: AbstractDataFrame, V <: FMExpr}

  #special case of no covariates
  if isnothing(rhs) #(rhs == Symbol(""))
    return ones(Float64,size(df,1),1)
  end

  #WARNING: This is a HACK! Replace with `nothing` or something else (See StatsModels github)
  #or just redo the whole thing in the GLM framework
  #Will be able to replace the LHS w/0 in subsequent versions
  lhs_hack::Symbol = names(df)[1]
  f::FormulaTerm = @eval(@formula($lhs_hack ~ $rhs))

  return generateX(M, df, f)
end

#helper function to create a one-sided formula given an expression
#IN: an expression and dataframe
#OUT: A one-sided formula object
#WARNING: Broken! if everything works drop this
#=function get1SidedFormula(RHS::T)::FormulaTerm where
  {T <: FMExpr}
  return @eval(@formula(nothing ~ $RHS)) #WARNING: THIS IS A HACK!!!!
end=#

#helper function  to get the list of symbols in an expression
#IN: A dataframe and a string object, typically representing a formula
#OUT: a vector of symbols in the string object
function getSymbolsFromExpr(inexpr::V)::Vector{Symbol} where {V <: FMExpr}
  isnothing(inexpr) && return Vector{Symbol}() #special case with no focal var
  symstrings::Vector{String} = split(string(inexpr),
    ['|','=','~',' ','+','*','&',')','(','-'], keepempty=false)

  filter!(s::String->typeof(Meta.parse(s))≠Int,symstrings)

  #convert the strings to symbols
  return Symbol.(symstrings)
end




#drops missings form a dataframe
function cloakmissings(::Type{M}, df::T, syms::Vector{Symbol})::SubDataFrame where {
    T<:AbstractDataFrame, M<:AbstractMatrix}
  local sdf::SubDataFrame

  (M <: CuArray) && CuArrays.allowscalar(true) #WARNING- this is expensive!
  sdf = view(df,completecases(df[!,syms]),:)
  (M <: CuArray) && CuArrays.allowscalar(false)

  return sdf
end

#helper function for the above, de-nulls entire dataframe
#IN: a subdataframe #OUT: subdataframe less the null entries
function cloakmissings(::Type{M}, df::T)::SubDataFrame where {T<:AbstractDataFrame, M<:AbstractMatrix}
  return cloakmissings(M, df, names(df))
end

#####################FMQR Type###########################
#=This is a storage object specific to the QR optimizations
meant to hold the QR factorization and a handle for X
Values: X matrix that is decomposed, Q matrix, and R matrix
The X matrix has dimensions of NxK. R is the inverse, stored
for computational efficiency.=#

struct FMQR{M<:AbstractMatrix}
  Q::M
  R::M

  N::Int
  K::Int
  Rinv::M
end

#main constructor for the FMQR object
#IN: X Matrix
#OUT: FMQR object
function FMQR(::Type{T}, X::M)::FMQR where {T<:AbstractMatrix, M<:AbstractMatrix}
  local Q::M
  local R::M
  local Rinv::M
  local QRcompact

  if M == T
    QRcompact = qr(X)
  else
    QRcompact = qr(T(X))
  end


  R = QRcompact.R

  #hack to get around a bug involving scalar operations when pulling the Q matrix
  Q = T<:CuArray ? CuMatrix(QRcompact.Q) : QRcompact.Q

  try
    Rinv = T<:CuArray ? Matrix(R)\I : R\I
  catch
    println("X matrix:")
    display(X)

    println("R matrix:")
    display(R)
    error("Inversion of R matrix failed!")
  end
  return FMQR{M}(Q,
      R,
      size(X,1),
      size(X,2),
      Rinv)
end

#default constructor, enables legacy compatability
function FMQR(X::M)::FMQR{M} where {M<:AbstractMatrix}
  return FMQR(M, X)
end


#####################FMLM Model Type#########################
#This is designed to hold the minimal info for a linear model
#Plus the associated FMQR object
struct FMLM{M<:AbstractMatrix, V<:AbstractVector} <: FMModel
  xqr::FMQR #handle to the xqr object
  X::M #the data and x matrix (handle)
  Y::V #Y data

  N::Int #Number of data points (rows)
  K::Int #Number of data columns
  dof::Int
  β::V #beta coefficient
  ε::V #residuals

  Xnames::Vector{Symbol} #names of the X variables
  Yname::Symbol #names of the Y variables

  clusters::Vector{C} where C <: FMData #holds the clusters

end


function FMLM(X::M, Y::V; clusters::Vector{C}=[nothing],
        Xnames::Vector{Symbol} = Symbol.(:X, 1:size(X,2)),
        Yname::Symbol = :Y, dof = size(X,1)-size(X,2),
        qrtype::Type = M)::FMLM where {
          M<:AbstractMatrix, V<:AbstractVector, C<:FMData}
    local xqr::FMQR

    xqr = FMQR(qrtype, X)

    β::V = Vector{Float64}(undef, xqr.K)
    coeff!(xqr,Y,β)
    ε::V = Vector{Float64}(undef, xqr.N)
    resid!(X,Y,β,ε)

    return FMLM{M,V}(xqr, X, Y, xqr.N, xqr.K, dof, β, ε, Xnames, Yname, clusters)
end


#=Constructor and helper method which gets required info from DataFrame
IN: The source dataframe, a formula expression for the RHS,
    the dependent variable as a symbol, optionally names for X columns,
    Y and a switch to elliminate null values
OUT: A FMLM object
NOTE: RHS expression must be ordered with factors last, otherwise
    the factor names will be incorrect =#
function FMLM(df::AbstractDataFrame,  Xexpr::FMExpr, Ysym::Symbol,
    ::Type{M} = Matrix{Float64}, ::Type{V} = Vector{Float64};
    Xnames::Vector{Symbol}=[Symbol("Intercept")], Yname::Symbol=:Y,
    containsmissings::Bool=true, withinsym::NSymbol = nothing, clustersym::NSymbol = nothing,
        qrtype::Type = M, checkwithin::Bool = false)::FMLM  where {
         M<:AbstractMatrix, V<:AbstractVector}

    local Xmodelmatrix::M
    local Y::V


    #get the list of symbols from the expression for X
    Xsyms::Vector{Symbol} = getSymbolsFromExpr(Xexpr)

    #println("flag2.2")
    #form an array of all symbols for which we want to drop nulls
    sarray::Vector{Symbol} = [Xsyms; Ysym]
    if withinsym ≠ nothing && withinsym ∉ sarray
      sarray = [sarray; withinsym]
    end
    if clustersym ≠ nothing && clustersym ∉ sarray
      sarray = [sarray; clustersym]
    end
    sdf::SubDataFrame = view(df[!, sarray], 1:size(df,1), :)
    #println("flag1")
    if containsmissings    # if we are going to drop nulls
      sdf = cloakmissings(M, sdf)
    end

    #get the within and clustering symbols
    if withinsym ≠ nothing
      within::Vector{W where W <: FMData}  = Vector{W where W <: FMData}(sdf[!, withinsym])
    end

    #println("flag4.2")
    if clustersym ≠ nothing
      clusters::Vector{C where C <: FMData}  = Vector{C where C <: FMData}(sdf[!, clustersym])
    else
      clusters = [nothing]
    end

    #println("flag2")
    #drop the columns which arn't going into the model matrix
    sdf = view(sdf, 1:size(sdf,1), [Xsyms; Ysym])
    #println("flag3")

    #get the factor expanded model matrix (Adds the dummies)
    Xmodelmatrix = generateX(M, sdf, Xexpr)
    Y = sdf[!, Ysym]

    #println("flag41")

    #this is a check to make sure the factor expansion follows other columns
    if length(Xsyms) > 1
        for i ∈ 2:length(Xsyms) #all types
            if typeof(sdf[!, Xsyms[i-1]])<:CategoricalArray &&
                !(typeof(sdf[!, Xsyms[i]])<:CategoricalArray)
                @warn "Type of $(typeof(sdf[i])) is preceeded by type of
                    factor. Column names are likely MISALIGNED. Put factors
                    after numerical variables in RHS expressions."
            end
        end
    end

    #Now assign the names as appropriate
    Xnamesfull::Vector{Symbol} = Vector{Symbol}(undef, size(Xmodelmatrix,2))
    for i ∈ 1:length(Xnamesfull)
        Xnamesfull[i] = i≤length(Xnames) ? Xnames[i] : Symbol(:X,i)
    end

    #get the appropriate matrices and construct the FMLM object

    if withinsym != nothing
        return FMLM(Xmodelmatrix, Y,
          within, clusters = clusters,
            Xnames=Xnamesfull, Yname=Yname, qrtype=qrtype)
    elseif (!checkwithin)
      #  println("flag3b")
        return FMLM(Xmodelmatrix,
            Y, clusters = clusters,
            Xnames=Xnamesfull, Yname=Yname, qrtype=qrtype)
    else
      #  println("flag3b")
        return FMLMold(Xmodelmatrix,
            Y, clusters = clusters,
            Xnames=Xnamesfull, Yname=Yname, qrtype=qrtype)
    end

end

#transforms the data using the within estimator, including an adjustment to dof
function FMLM(X::M, Y::V, within::Vector{W};
    clusters::Vector{C} = [nothing],
    Xnames::Vector{Symbol} = (Symbol).(:X, 1:size(X,2)),
    Yname::Symbol = :Y,
    qrtype::Type = M)::FMLM where {
      W<:FMData, C<:FMData, M<:AbstractMatrix, V<:AbstractVector}

    #println("flag4")
    #we get the unique values and make a two way table
    #tVars::Vector{T} = unique(tFI)
    ivars::Vector{W} = unique(within)
    iN::Int = length(ivars)
    itable::Dict{W,Int} = Dict{W,Int}(ivars[i] => i for i::Int ∈ 1:iN)
    iindex::Dict{W,Vector{Bool}} = Dict(iv=> within.==iv for iv ∈ ivars)
    iXindex::Dict{W,SubArray} = Dict(iv=>view(X, iindex[iv], :) for iv ∈ ivars)
    iYindex::Dict{W,SubArray} = Dict(iv=>view(Y, iindex[iv]) for iv ∈ ivars)

    K::Int = size(X,2)
    N::Int = size(X,1)

    #now create the codified versions
    #Xi::Matrix{Float64} = Matrix{Float64}(length(ivars),size(X,2))
    withincode::Vector{Int} = ((x::W)->itable[x]).(within)

    meanXY::Matrix{Float64} = zeros(iN,K+1)
    NXY::Vector{Int} = zeros(iN)

    #println("$N,$K,$iN,$(size(NXY))")
    @fastmath for (i::Int, iv::W) ∈ enumerate(ivars)
      meanY::Float64 = mean(iYindex[iv])
      meanX::M = mean(iXindex[iv], dims=1) #This is actually a row vector
      meanX[1] = 0.0 #WARNING: Assumes first column is the intercept

      iYindex[iv] .-= meanY
      iXindex[iv] .-= meanX
    end


    return FMLM(X, Y,
      clusters=clusters, Xnames=Xnames, Yname=Yname,
      dof = N - iN - K, qrtype=qrtype)
end

function FMLMold(X::M, Y::V, within::AbstractVector{W};
    clusters::Vector{C} = [nothing],
    Xnames::Vector{Symbol} = (Symbol).(:X, 1:size(X,2)),
    Yname::Symbol = :Y,
    qrtype::Type = M)::FMLM where {
      W<:FMData, C<:FMData, M<:AbstractMatrix, V<:AbstractVector}

    if M<:CuArray
      error("Cluster code needs overhaul for effective gpu usage")
    end


    #println("flag4")
    #we get the unique values and make a two way table
    #tVars::Vector{T} = unique(tFI)
    ivars::Vector{W} = unique(within)
    iN::Int = length(ivars)
    itable::Dict{W,Int} = Dict{W,Int}(ivars[i] => i for i::Int ∈ 1:iN)

    K::Int = size(X,2)
    N::Int = size(X,1)

    #now create the codified versions
    #Xi::Matrix{Float64} = Matrix{Float64}(length(ivars),size(X,2))
    withincode::Vector{Int} = ((x::W)->itable[x]).(within)

    meanXY::Matrix{Float64} = zeros(iN,K+1)
    NXY::Vector{Int} = zeros(iN)

    #println("$N,$K,$iN,$(size(NXY))")
    @fastmath for r::Int ∈ 1:N
      meanXY[withincode[r],K+1] += Y[r]
      NXY[withincode[r]] += 1.0
    end

    @fastmath for c::Int ∈2:K, r ∈1:N
      meanXY[withincode[r],c] += X[r,c]
    end

    meanXY ./= NXY

    @fastmath for c::Int ∈2:K, r::Int ∈1:N
      X[r,c] -= meanXY[withincode[r],c]
    end

    @fastmath for r::Int ∈ 1:N
      Y[r] -= meanXY[withincode[r],K+1]
    end

    return FMLM(X, Y,
      clusters=clusters, Xnames=Xnames, Yname=Yname,
      dof = N - iN - K, qrtype=qrtype)
end

#this function calculates the Pearson correlation coefficient of the predicted values
#relative to the realized values
getR(args...; keyargs...) = error("use R instead")
R(lin::FMLM)::Float64 = cor(lin.Y, lin.Y .- lin.ε)

getR²(args...; keyargs...) = error("use R² instead")
function R²(lin::FMLM; adjusted::Bool = true)::Float64
  R²::Float64 = R(lin)^.2
  R² = adjusted ? R² - (1-R²)*(lin.N-lin.dof+1)/lin.dof : R²

  return R²
end


#fucntion for getting regression coefficients
getTerm(args...; keyargs...) = error("use term instead")
term(lm::FMLM, s::Symbol) = (lm.β[findfirst(lm.Xnames, s)])::Float64


#########################Project!##############################
#=gets the projection matrix

IN: THe FMQR decomposition and the memory for  the projection matrix
OUT: Writes and returns the projection matrix
NOTE: This method will allocate the projection matrix if no second
argument is supplied. =#
function project!(xqr::FMQR{M}, P::M=M(undef, xqr.N, xqr.N)
  )::Nothing where M<:AbstractMatrix

  T::Type = eltype(M)
  BLAS.gemm!('N','T',T(1.0),xqr.Q,xqr.Q,T(0.0),P)
  #use the BLAS library for fast matrix algebra
  return nothing

end

#=This version gets only the diagonal of the projection matrix
IN: THe FMQR decomposition and the memory for  the projection diagonal
OUT: Writes and returns the projection matrix=#
function project!(xqr::FMQR{M}, P::V)::Nothing where {M<:AbstractMatrix, V<:AbstractVector}
  #get the diagonal quickly
  #equivlenet to diag(Q*Q')
  P .= diag(xqr.Q*xqr.Q')

  return nothing

end

#=This version runs either of the above methods after generating
the xqr object given an independent variable matrix
IN: Independent variable X matrix and the memory for  the projection matrix
or vector for the diagonal
OUT: Writes and returns the projection matrix=#
function project!(X::M,
  P::AbstractArray=M(undef, size(X,1),size(X,1)))::Nothing where M<:AbstractMatrix

  return project!(FMQR(X),P)
end

#=simple test function to verify the projection optimizations
output should be identical to the matrix function above
IN: Independent variable X matrix, memory for the projection matrix
OUT: Writes and returns teh projection matrix=#
projectSlow!(args...; keyargs...) = error("Use projectslow! instead")
function projectslow!(X::M, P::AbstractArray) where M<:AbstractMatrix
  XXinv = M(Matrix(X' * X)\I)
  P .= X * XXinv * X'
  return P
end

#########################resid!#########################


#=Main method to get residuals
IN: X matrix the RHS variables, LHS dependent variable,
a coefficient vector, optionally the memory for the residual vector
OUT: Writes and returns the residual vector=#
getResid!(args...; keyargs...) = error("use resid! instead")
function resid!(X::M, Y::V,
  β::V, ε::V=similar(Y)) where {M<:AbstractMatrix, V<:AbstractVector}

  T::Type = eltype(V)
  return BLAS.axpy!(T(1.0), Y, BLAS.gemv!('N',T(-1.0), X, β,T(0.0),ε)) #ε=Y-Xβ
end

#########################getCoef!###########################
# Gets the coefficients in a regression



#=Main method to get the regression coefficients
IN: QR decomposition of X matrix the RHS variables, LHS dependent variable,
optionally the memory for the coefficient vector
OUT: Writes and returns the coefficient vector=#

getCoeff!(args...; keyargs...) = error("Use coeff! instead")

function coeff!(xqr::FMQR{M}, Y::V,
  β::V = V(undef, xqr.K))::V where  {M<:AbstractMatrix, V<:AbstractVector}

  #println("got here 5")
  #println("xqr.Rinv size: $(size(xqr.Rinv)) xqr.Q size:$(size(xqr.Q))")
  RinvQt::M = BLAS.gemm('N', 'T', xqr.Rinv, xqr.Q)
  #println("got here 6")

  if M<:Matrix{Float64}
    return  BLAS.gemv!('N',1.0,RinvQt,Y,0.0,β)
  else
    T::Type = eltype(M)
    return  BLAS.gemv!('N',T(1.0),RinvQt,Y,T(0.0),β)
  end
end


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

  Σ .*= dofcorrect #WARNING: dof correction? is it being double counted?

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
function homoskedasticΣslow(X::M, ε::V)::M where {M<:AbstractMatrix, V<:AbstractVector}
  T::Type = eltype(V)

  out = Matrix(X' * X)\I .* (ε⋅ε)/(size(X,1)-size(X,2))
  return out
end

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
homoskedasticΣslow(lin::FMLM{M where M<:AbstractMatrix}) =
  homoskedasticΣslow(lin.X, lin.Y.-lin.X * lin.β)::M


###################getWhiteERRORs ###OLS only

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

  return XXinv * X' * diagm(ε) .^ 2.0 * X * XXinv
end

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
whiteΣslow(lin::FMLM) = whiteΣslow(lin.X, lin.Y.-lin.X * lin.β, lin.N/lin.dof)

###########################getClustered SEs
#method for getting clustered errors
#from http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf
#forumala: [X'X]^-1*B*[X'X]^-1 where B= sum over G (Xg'*εg*εg'*Xg) and g is indexed for G clusters
getClustered!(args...; keyargs...) = error("use clusteredΣ! instead")
function clusteredΣ!(X::M, xqr::FMQR{M}, ε::V, clusters::C,
    Σ::M = M(undef, xqr.K, xqr.K),
    dof::Int = xqr.N-xqr.K)::M where {M<:AbstractMatrix, V<:AbstractVector, C<:AbstractVector}

  #convenience vectors and values
  cvars::C = unique(clusters)
  Nclusters::Int = length(cvars)
  ctable::Dict = Dict(cvars[i] => i for i::Int ∈ 1:Nclusters)
  clustercode::Vector{Int} = ((x)->ctable[x]).(clusters)

  #iterate through the groups
  B::M = zeros(xqr.K, xqr.K)
  @fastmath for i::Int ∈ 1:Nclusters
    Xg::SubArray = view(X, clustercode .== i, :)
    ug::SubArray = view(ε, clustercode .== i)
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

  #=@fastmath for i ∈ 1:xqr.K
    if Σ[i,i] < 0.
      error("Negative variance ($i): $(diag(Σ))")
    end
  end=#

  #println("$(diag(Σ))")
  return Σ
end

function clusteredΣ!(lin::FMLM{M}, Σ::M = M(undef, lin.K, lin.K))::M where M<:AbstractMatrix

  return clusteredΣ!(lin.X, lin.xqr, lin.ε, lin.clusters, Σ, lin.dof)
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


  function neweywestΣ!(lin::FMLM, lag::Int,
    Σ::M = M(undef, lin.K, lin.K))::M where M<: AbstractMatrix

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
function neweywestpanelΣ!(X::M, xqr::FMQR{M}, ε::V, clusters::Vector{C},
    lag::Int, Σ::M = M(undef, xqr.K, xqr.K))::M where {
    M<:AbstractMatrix, V<:AbstractVector, C}

  N::Int = length(ε)
  E::Type = eltype(M)

  K::Int = size(X,2)
  Sₜ::M = zeros(K,K) #holds teh central matrix

  temp::M = M(undef, K, K) #pre-allocate working matrix
  RRinv::M = BLAS.gemm('N', 'T', xqr.Rinv, xqr.Rinv) #this is equivelent to [X'X]^-1

  cvars::Vector{C} = unique(clusters)
  cindex::Dict = Dict(cvar => clusters .== cvar for cvar ∈ cvars)

  #iterate over all stocks
  for cvar ∈ cvars
    Xₙ::SubArray = view(X, cindex[cvar], :)
    εₙ::SubArray = view(ε, cindex[cvar])

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


neweywestΣfunc(lag::Int) = (lin::FMLM)-> neweywestpanelΣ!(lin, lag)


function neweywestpanelΣslow!(X::M, xqr::FMQR{M}, ε::V, clusters::Vector{C},
    lag::Int, Σ::M = M(undef, xqr.K, xqr.K))::M where {
    M<:AbstractMatrix, V<:AbstractVector, C}

  N::Int = size(X,1)
  K::Int = size(X,2)
  T::Type = eltype(M)

  #pre-allocate
  Sₜ::M = zeros(K, K)

  #this method requires a sorted array
  #will do this via dataframes
  df::DataFrame = DataFrame(deepcopy(X))
  xnames = names(df)
  df.clusters = deepcopy(clusters)
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
