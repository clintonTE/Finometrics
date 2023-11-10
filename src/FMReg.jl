
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


  m::M = ModelMatrix(ModelFrame(f, df)).m

  return m
end

#Same as above but allows for an expression
function generateX(::Type{M}, df::T, rhs::V)::M where
  {M<:AbstractMatrix, T <: AbstractDataFrame, V <: FMExpr}

  #special case of no covariates
  if isnothing(rhs) #(rhs == Symbol(""))
    return ones(Float64,size(df,1),1)
  end

  f::FormulaTerm = @eval(@formula(0 ~ $rhs))

  return generateX(M, df, f)
end

#helper function  to get the list of symbols in an expression
#IN: A dataframe and a string object, typically representing a formula
#OUT: a vector of symbols in the string object
function symbolsFromExpr(inexpr::V)::Vector{Symbol} where {V <: FMExpr}
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

  sdf = view(df,completecases(df[!,syms]),:)

  return sdf
end

#helper function for the above, de-nulls entire dataframe
#IN: a subdataframe #OUT: subdataframe less the null entries
function cloakmissings(::Type{M}, df::T)::SubDataFrame where {T<:AbstractDataFrame, M<:AbstractMatrix}
  return cloakmissings(M, df, propertynames(df))
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
  #hack disabled
  #Q = T<:CuArray ? CuMatrix(QRcompact.Q) : QRcompact.Q
  Q = QRcompact.Q
  

  try
    #hack below disabled
    #Rinv = T<:CuArray ? pinv(Matrix(R)) : pinv(R)
    Rinv = pinv(R)

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

const FMClusters = Vector{Vector{C} where C<:FMData}

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

  Xnames::Vector{String} #names of the X variables
  Yname::String #names of the Y variables

  clusters::FMClusters  #holds the clusters

end


function FMLM(X::M, Y::V; clusters::FMClusters = FMClusters([]),
        Xnames::Vector{String} = string.(:X, 1:size(X,2)),
        Yname::String = "Y", dof = size(X,1)-size(X,2),
        qrtype::Type = M)::FMLM where {
          M<:AbstractMatrix, V<:AbstractVector}
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
    Xnames::Vector{String}=Vector{String}(), Yname::String="Y",
    containsmissings::Bool=true, withinsym::NSymbol = nothing,
    clustersyms::C= nothing,
    qrtype::Type = M, checkwithin::Bool = false)::FMLM  where {
     M<:AbstractMatrix, V<:AbstractVector, C<:Union{NSymbol, Vector{Symbol}}}

    local Xmodelmatrix::M
    local Y::V
    local clusters::FMClusters = []

    (C <: Vector) && (length(clustersyms)>2) && error("Max of 2 clusters supported!!!")

    #get the list of symbols from the expression for X
    Xsyms::Vector{Symbol} = symbolsFromExpr(Xexpr)

    #println("flag2.2")
    #form an array of all symbols for which we want to drop nulls
    sarray::Vector{Symbol} = [Xsyms; Ysym]
    if withinsym ≠ nothing && withinsym ∉ sarray
      sarray = [sarray; withinsym]
    end
    if !isnothing(clustersyms)
      for clustersym ∈ [clustersyms;]
        (clustersym ∉ sarray) && (sarray = [sarray; clustersym])
      end
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
    if !isnothing(clustersyms)
      for clustersym ∈ [clustersyms;]
        push!(clusters, Vector(sdf[!, clustersym]))
      end
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
    Xnamesfull::Vector{String} = Vector{String}(undef, size(Xmodelmatrix,2))
    for i ∈ 1:length(Xnamesfull)
        Xnamesfull[i] = i≤length(Xnames) ? Xnames[i] : string(:X,i)
    end

    #get the appropriate matrices and construct the FMLM object

    if (withinsym != nothing) && (!checkwithin)
        return FMLM(Xmodelmatrix, Y,
          within, clusters = clusters,
            Xnames=Xnamesfull, Yname=Yname, qrtype=qrtype)
    elseif (withinsym != nothing) && (checkwithin)
      #  println("flag3b")
        return FMLMaltwithin(Xmodelmatrix,
            Y, within, clusters = clusters,
            Xnames=Xnamesfull, Yname=Yname, qrtype=qrtype)
    else
      #println("flag3b")
      return FMLM(Xmodelmatrix,
          Y, clusters = clusters,
          Xnames=Xnamesfull, Yname=Yname, qrtype=qrtype)
    end

end


#transforms the data using the within estimator, including an adjustment to dof
function FMLM(X::M, Y::V, within::AbstractVector{W};
    clusters::FMClusters = FMClusters([]),
    Xnames::Vector{String} = (string).(:X, 1:size(X,2)),
    Yname::String = "Y",
    qrtype::Type = M)::FMLM where {
      W<:FMData, M<:AbstractMatrix, V<:AbstractVector}

    Xmat::Matrix{Float64} = X
    Yvec::Vector{Float64} = Y
    Xstart::Int = 1

    if minimum(Xmat[:,1]) ≈ maximum(Xmat[:,1]) ≈ 1.0
      @warn("intercept with within effects detected- problem is not well posed
        answer may be wrong!!! (it might be ok, but why take the risk?)
        remediate by using +0 in the rhs regression expression")
      Xstart = 2
    end

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
      meanXY[withincode[r],K+1] += Yvec[r]
      NXY[withincode[r]] += 1.0
    end

    @fastmath for c::Int ∈Xstart:K, r ∈1:N
      meanXY[withincode[r],c] += Xmat[r,c]
    end

    meanXY ./= NXY

    @fastmath for c::Int ∈Xstart:K, r::Int ∈ 1:N
      Xmat[r,c] -= meanXY[withincode[r],c]
    end

    @fastmath for r::Int ∈ 1:N
      Yvec[r] -= meanXY[withincode[r],K+1]
    end

    X = Xmat
    Y = Yvec

    return FMLM(X, Y,
      clusters=clusters, Xnames=Xnames, Yname=Yname,
      dof = N - iN - K, qrtype=qrtype)
end

#this alternate method is easier to understand, but as of 1.2 it is 20% slower
#I expect the speed of this version to increase over time
function FMLMaltwithin(X::M, Y::V, within::Vector{W};
    clusters::FMClusters = FMClusters([]),
    Xnames::Vector{String} = (string).(:X, 1:size(X,2)),
    Yname::String = "Y",
    qrtype::Type = M)::FMLM where {
      W<:FMData, M<:AbstractMatrix, V<:AbstractVector}

    #println("flag4")
    #we get the unique values and make a two way table
    #tVars::Vector{T} = unique(tFI)
    ivars::Vector{W} = unique(within)
    iN::Int = length(ivars)

    #perposefully do not use a copy command, since we are looking to replace original elements
    Xmat::Matrix{Float64} = X
    Yvec::Vector{Float64} = Y

    K::Int = size(X,2)
    N::Int = size(X,1)

    firstcolintercept::Bool = false
    if minimum(Xmat[:,1]) ≈ maximum(Xmat[:,1]) ≈ 1.0
      @warn("intercept with within effects detected- problem is not well posed
        answer may be wrong!!! (it might be ok, but why take the risk?)
        remediate by using +0 in the rhs regression expression")
      firstcolintercept = true
    end


    #println("$N,$K,$iN,$(size(NXY))")
    #warning change these if multi-threading
    meanX::Matrix{Float64} = Matrix{Float64}(undef, 1,K) #WARNING:
    local sX::SubArray
    local sY::SubArray
    iindex::Vector{Bool} = Vector{Bool}(undef, length(within))

    @fastmath for (i::Int, iv::W) ∈ enumerate(ivars)

      iindex .= within.==ivars[i]
      sX=  view(Xmat, iindex, :)
      sY =  view(Yvec, iindex)

      meanY::Float64 = mean(sY)
      meanX .= mean(sX, dims=1) #This is actually a row vector

      if firstcolintercept
        meanX[1] = 0.0 #WARNING: Assumes first column is the intercept
      end

      sY .-= meanY
      sX .-= meanX
    end

    X = Xmat
    Y = Yvec


    return FMLM(X, Y,
      clusters=clusters, Xnames=Xnames, Yname=Yname,
      dof = N - iN - K, qrtype=qrtype)
end

#this function calculates the Pearson correlation coefficient of the predicted values
#relative to the realized values
getR(args...; keyargs...) = error("use R instead")
R(lin::FMLM)::Float64 = cor(lin.Y, lin.Y .- lin.ε)

getR²(args...; keyargs...) = error("use R² instead")

function R²(lin::FMLM; adjusted::Bool = false)::Float64
  R²::Float64 = R(lin)^2
  R² = adjusted ? R² - (1-R²)*(lin.N-lin.dof+1)/lin.dof : R²

  return R²
end


#fucntion for getting regression coefficients
getTerm(args...; keyargs...) = error("use βterm instead")
βterm(lm::FMLM, s::String) = (lm.β[findfirst(lm.Xnames, s)])::Float64
βterm(lm::FMLM, s::Symbol) = βterm(lm, string(s))


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
  P .= vec(sum(xqr.Q .^ 2, dims=2))

  return nothing

end

#=This version runs either of the above methods after generating
the xqr object given an independent variable matrix
IN: Independent variable X matrix and the memory for  the projection matrix
or vector for the diagonal
OUT: Writes and returns the projection matrix=#
function project!(X::M,
  P::AbstractArray=M(undef, size(X,1),size(X,1)))::Nothing where M<:AbstractMatrix

  return project!(FMQR(M, X),P)
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
