
##################Utility Functions ###############################

#helper function  to get the model matrix from a dataframe given a formula
#IN: A dataframe and formula object
#OUT: the model matrix
function getModelMatrix(df::T, f::FormulaTerm)::Matrix{Float64} where
  T <: AbstractDataFrame
  return ModelMatrix(ModelFrame(f, df, model=LinearModel)).m
end

#Same as above but allows for an expression
function getModelMatrix(df::T, exp::V)::Matrix{Float64} where
  {T <: AbstractDataFrame, V <: FMExpr}

  #special case which crashes Formula
  if exp == Symbol("")
    return ones(Float64,size(df,1),1)
  end

  return getModelMatrix(df, get1SidedFormula(exp))
end


#helper function to create a one-sided formula given an expression
#IN: an expression and dataframe
#OUT: A one-sided formula object
function get1SidedFormula(RHS::T)::FormulaTerm where
  {T <: FMExpr}
  return @eval(@formula(identity(1) ~ $RHS)) #WARNING: THIS IS A HACK!!!!
end

#helper function  to get the list of symbols in an expression
#IN: A dataframe and a string object, typically representing a formula
#OUT: a vector of symbols in the string object
function getSymbolsFromExpr(exp::V)::Vector{Symbol} where {V <: FMExpr}
  symStrings::Vector{String} = split(string(exp),
    ['|','=','~',' ','+','*','&',')','(','-'], keepempty=false)

  filter!(s::String->typeof(Meta.parse(s))≠Int,symStrings)

  #convert the strings to symbols
  return Symbol.(symStrings)
end




#drops missings form a dataframe
function cloakmissings(df::T, syms::Vector{Symbol})::T where {T<:AbstractDataFrame}
  return view(df,completecases(df[!,syms]),:)
end

#helper function for the above, de-nulls entire dataframe
#IN: a subdataframe #OUT: subdataframe less the null entries
function cloakmissings(df::T)::T where {T<:AbstractDataFrame}
  return cloakmissings(df, names(df))
end

#####################FMQR Type###########################
#=This is a storage object specific to the QR optimizations
meant to hold the QR factorization and a handle for X
Values: X matrix that is decomposed, Q matrix, and R matrix
The X matrix has dimensions of NxK. R is the inverse, stored
for computational efficiency.=#

struct FMQR
  X::Matrix{Float64}
  Q::Matrix{Float64}
  R::Matrix{Float64}

  N::Int
  K::Int
  RInv::Matrix{Float64}
end

#helper method which takes the inverse of the R matrix (which is used in a lot)
#IN: X matrix and its Q R decomposition in matrix form
#OUT: FMQR Object
function FMQR(X::Matrix{Float64},Q::Matrix{Float64},R::Matrix{Float64};keepX=true)::FMQR
  #=println("Keep: $keepX")
  println("X Size: $(size(X))")
  println("Q Size: $(size(Q))")
  println("R Size: $(size(R))")=#
  #display(
  FMQR(keepX ? X : Matrix{Float64}(),
    Q,
    R,
    size(X,1),
    size(X,2),
    (R)\Matrix{Float64}(I,size(X,2),size(X,2)))
end

#main constructor for the FMQR object
#IN: X Matrix
#OUT: FMQR object
function FMQR(X::Matrix{Float64}; keepX::Bool = true)::FMQR
  local Q::Matrix{Float64}
  local R::Matrix{Float64}
  local QRCompact::LinearAlgebra.QRCompactWY

  try
    QRCompact = qr(X)
  catch err
    display(X)
    error("Error with qr decomposition msg2389 error: $err")
  end

  Q = Matrix(QRCompact.Q)
  R = Matrix(QRCompact.R)

  return FMQR(X,Q,R)
end

#default constructor to clear memory
FMQR() = FMQR(Matrix{Float64}(),Matrix{Float64}(),Matrix{Float64}(),
    0,0,Matrix{Float64}())::FMQR


#####################FMLM Model Type#########################
#This is designed to hold the minimal info for a linear model
#Plus the associated FMQR object
 mutable struct FMLM <: FMModel
  xqr::FMQR #handle to the xqr object
  X::Matrix{Float64} #the data and x matrix (handle)
  Y::Vector{Float64} #Y data

  N::Int #Number of data points (rows)
  K::Int #Number of data columns
  dof::Int
  β::Vector{Float64} #beta coefficient
  ε::Vector{Float64} #residuals

  XNames::Vector{Symbol} #names of the X variables
  YName::Symbol #names of the Y variables

  clusters::Vector{C} where C <: FMData #holds the clusters

end


function FMLM(X::Matrix{Float64}, Y::Vector{Float64}; clusters::Vector{C}=[nothing],
        XNames::Vector{Symbol} = Symbol.(:X, 1:size(X,2)),
        YName::Symbol = :Y, dof = size(X,1)-size(X,2))::FMLM where C<:FMData

    xqr::FMQR = FMQR(X, keepX=false)
    β::Vector{Float64} = Vector{Float64}(undef, xqr.K)
    getCoeff!(xqr,Y,β)
    ε::Vector{Float64} = Vector{Float64}(undef, xqr.N)
    getResid!(X,Y,β,ε)

    return FMLM(xqr, X, Y, xqr.N, xqr.K, dof, β, ε, XNames, YName, clusters)
end

#=Constructor and helper method which gets required info from DataFrame
IN: The source dataframe, a formula expression for the RHS,
    the dependent variable as a symbol, optionally names for X columns,
    Y and a switch to elliminate null values
OUT: A FMLM object
NOTE: RHS expression must be ordered with factors last, otherwise
    the factor names will be incorrect =#
function FMLM(df::AbstractDataFrame,  XExpr::T, YSym::Symbol;
    XNames::Vector{Symbol}=[Symbol("Intercept")], YName::Symbol=:Y,
    containsmissings::Bool=true, withinSym::V = nothing, clusterSym::R = nothing,
        fixedEffectsSym::NSymbol = nothing)::FMLM  where {
        T <: FMExpr, V<:Union{Symbol,Nothing}, R<:Union{Symbol, Nothing}}

    #NOTE: Delete when code is fully updated
    if fixedEffectsSym ≠ nothing
      error("fixedEffectsSym no longer used. Replace within withinSym and/or clusterSym")
    end

    #get the list of symbols from the expression for X
    XSym::Vector{Symbol} = getSymbolsFromExpr(XExpr)

    #println("flag2.2")
    #form an array of all symbols for which we want to drop nulls
    viewArray::Vector{Symbol} = [XSym; YSym]
    if withinSym ≠ nothing && withinSym ∉ viewArray
      viewArray = [viewArray; withinSym]
    end
    if clusterSym ≠ nothing && clusterSym ∉ viewArray
      viewArray = [viewArray; clusterSym]
    end
    dfSub::SubDataFrame = view(df[!, viewArray], 1:size(df,1), :)
    #println("flag1")
    if containsmissings    # if we are going to drop nulls
      dfSub = cloakmissings(dfSub)
    end

    #get the within and clustering symbols
    if withinSym ≠ nothing
      within::Vector{W where W <: FMData}  = Vector{W where W <: FMData}(dfSub[withinSym])
    end

    #println("flag4.2")
    if clusterSym ≠ nothing
      clusters::Vector{C where C <: FMData}  = Vector{C where C <: FMData}(dfSub[clusterSym])
    else
      clusters = [nothing]
    end

    #println("flag2")
    #drop the columns which arn't going into the model matrix
    dfSub = view(dfSub, 1:size(dfSub,1), [XSym; YSym])
    #println("flag3")

    #get the factor expanded model matrix (Adds the dummies)
    XModelMatrix::Matrix{Float64} = getModelMatrix(dfSub, XExpr)
    #println("flag41")

    #this is a check to make sure the factor expansion follows other columns
    if length(XSym) > 1
        for i ∈ 2:length(XSym) #all types
            if typeof(dfSub[!, XSym[i-1]])<:CategoricalArray &&
                !(typeof(dfSub[!, XSym[i]])<:CategoricalArray)
                @warn "Type of $(typeof(dfSub[i])) is preceeded by type of
                    factor. Column names are likely MISALIGNED. Put factors
                    after numerical variables in RHS expressions."
            end
        end
    end

    #Now assign the names as appropriate
    XNamesFull::Vector{Symbol} = Vector{Symbol}(undef, size(XModelMatrix,2))
    for i ∈ 1:length(XNamesFull)
        XNamesFull[i] = i≤length(XNames) ? XNames[i] : Symbol(:X,i)
    end

    #get the appropriate matrices and construct the FMLM object

    if withinSym != nothing
      #println("flag3a")
            #println(typeof(Vector(dfSub[fixedEffectsSym])))
        return FMLM(XModelMatrix, Vector{Float64}(dfSub[YSym]),
          within, clusters = clusters,
            XNames=XNamesFull, YName=YName)
    else
      #  println("flag3b")
        return FMLM(XModelMatrix,
            Vector{Float64}(dfSub[!, YSym]), clusters = clusters,
            XNames=XNamesFull, YName=YName)
    end

end

#transforms the data using the within estimator, including an adjustment to dof
function FMLM(X::Matrix{Float64}, Y::Vector{Float64}, within::Vector{W};
    clusters::Vector{C} = [nothing],
    XNames::Vector{Symbol} = (Symbol).(:X, 1:size(X,2)),
    YName::Symbol = :Y)::FMLM where {W<:FMData, C<:FMData}


    #println("flag4")
    #we get the unique values and make a two way table
    #tVars::Vector{T} = unique(tFI)
    iVars::Vector{W} = unique(within)
    iN::Int = length(iVars)
    iTable::Dict = Dict(iVars[i] => i for i::Int ∈ 1:iN)

    K::Int = size(X,2)
    N::Int = size(X,1)

    #now create the codified versions
    #Xi::Matrix{Float64} = Matrix{Float64}(length(iVars),size(X,2))
    withinCode::Vector{Int} = ((x::W)->iTable[x]).(within)

    meanXY::Matrix{Float64} = zeros(iN,K+1)
    NXY::Vector{Int} = zeros(iN)

    #println("$N,$K,$iN,$(size(NXY))")
    @fastmath for r::Int ∈ 1:N
      meanXY[withinCode[r],K+1] += Y[r]
      NXY[withinCode[r]] += 1.0
    end

    @fastmath for c::Int ∈2:K, r ∈1:N
      meanXY[withinCode[r],c] += X[r,c]
    end

    meanXY ./= NXY

    @fastmath for c::Int ∈2:K, r::Int ∈1:N
      X[r,c] -= meanXY[withinCode[r],c]
    end

    @fastmath for r::Int ∈ 1:N
      Y[r] -= meanXY[withinCode[r],K+1]
    end

    return FMLM(X, Y, clusters=clusters, XNames=XNames, YName=YName, dof = N - iN - K)
end

#this function calculates the Pearson correlation coefficient of the predicted values
#relative to the realized values
getR(lin::FMLM)::Float64 = cor(lin.Y, lin.Y .- lin.ε)

function getR²(lin::FMLM; adjusted::Bool = true, percentage::Bool = false)::Float64
  R²::Float64 = getR(lin)^.2
  R² = adjusted ? R² - (1-R²)*(lin.N-lin.dof+1)/lin.dof : R²
  R² = percentage ? R²*100 : R²

  return R²
end


#fucntion for getting regression coefficients
getTerm(lm::FMLM, s::Symbol) = (lm.β[findfirst(lm.XNames, s)])::Float64


#########################Project!##############################
#=gets the projection matrix

IN: THe FMQR decomposition and the memory for  the projection matrix
OUT: Writes and returns the projection matrix
NOTE: This method will allocate the projection matrix if no second
argument is supplied. =#
function project!(xqr::FMQR, P::Matrix{Float64}=Matrix{Float64}(undef, xqr.N, xqr.N)
  )::Matrix{Float64}

  #use the BLAS library for fast matrix algebra
  return BLAS.gemm!('N','T',1.0,xqr.Q,xqr.Q,0.0,P)

end

#=This version gets only the diagonal of the projection matrix
IN: THe FMQR decomposition and the memory for  the projection diagonal
OUT: Writes and returns the projection matrix=#
function project!(xqr::FMQR, P::Vector{Float64})::Vector{Float64}
  #get the diagonal quickly
  #equivlenet to diag(Q*Q')
  P .= 0.0
  @fastmath for j::Int ∈ 1:xqr.K, i::Int ∈ 1:xqr.N
    P[i] += xqr.Q[i,j] * xqr.Q[i,j]
  end

  return P

end

#=This version runs either of the above methods after generating
the xqr object given an independent variable matrix
IN: Independent variable X matrix and the memory for  the projection matrix
or vector for the diagonal
OUT: Writes and returns the projection matrix=#
function project!(X::Matrix{Float64},
  P::Array{Float64}=Matrix{Float64}(undef, size(X,1),size(X,1)))::Array{Float64}

  return project!(FMQR(X),P)
end

#=simple test function to verify the projection optimizations
output should be identical to the matrix function above
IN: Independent variable X matrix, memory for the projection matrix
OUT: Writes and returns teh projection matrix=#
function projectSlow!(X::Matrix{Float64},P::Matrix{Float64})
  P .= (X * ((X' * X)\Matrix{Float64}(I,size(X,2),size(X,2))) * X')
  return P
end

#########################getResid!#########################


#=Main method to get residuals
IN: X matrix the RHS variables, LHS dependent variable,
a coefficient vector, optionally the memory for the residual vector
OUT: Writes and returns the residual vector=#
function getResid!(X::Matrix{Float64}, Y::T,
  β::T, ε::T=similar(Y))::T where T<:StridedVector{Float64}

  return BLAS.axpy!(1.0, Y, BLAS.gemv!('N',-1.0, X, β,0.0,ε)) #ε=Y-Xβ
end

#########################getCoef!###########################
# Gets the coefficients in a regression



#=Main method to get the regression coefficients
IN: QR decomposition of X matrix the RHS variables, LHS dependent variable,
optionally the memory for the coefficient vector
OUT: Writes and returns the coefficient vector=#
function getCoeff!(xqr::FMQR, Y::T,
  β::T = Vector{Float64}(undef, xqr.K))::T where  T<:StridedVector{Float64}

  #println("got here 5")
  #println("xqr.RInv size: $(size(xqr.RInv)) xqr.Q size:$(size(xqr.Q))")
  M::Matrix{Float64} = BLAS.gemm('N', 'T', xqr.RInv, xqr.Q)
  #println("got here 6")

  return  BLAS.gemv!('N',1.0,M,Y,0.0,β)
end

#=Convenience method to get the regression coefficients
IN: X matrix the RHS variables, LHS dependent variable,
optionally the memory for the coefficient vector
OUT: Writes and returns the coefficeint vector=#
getCoeff!(X::Vector{Float64}, Y::Vector{Float64}) =
  getCoeff!(FMQR(X), Y, Vector{Float64}(undef, size(X,2)))::Vector{Float64}


###################getHomosked!#######################
#Gets the covariance matrix under Homoskedastic assumptions σ^2[X'X]^-1

#main method to get the covariance matrix
#IN: The QR decomposition from the regression, a vector of residuals,
#optionally memory for the covariance matrix
#OUT: Writes and returns the covariance matrix
function getHomoskedΣ!(xqr::FMQR, ε::Vector{Float64},
  Σ::Matrix{Float64} = Matrix{Float64}(undef, lin.K, lin.K); dofCorrect::Float64 = 1.0)::Matrix{Float64}

  return dofCorrect .* BLAS.gemm!('N','T',(ε⋅ε)/(xqr.N-xqr.K),xqr.RInv,xqr.RInv,0.0,Σ) #[X'X]^-1*σ2
end

#=helper method for the above function which extracts the QR decomposition
from a linear model
IN: A linear model, optionally memory for the covariance matrix
OUT: WRites and returns the covariance matrix=#
function getHomoskedΣ!(lin::FMLM,
  Σ::Matrix{Float64}=Matrix{Float64}(undef, lin.K, lin.K))::Matrix{Float64}

  return getHomoskedΣ!(lin.xqr, lin.ε, Σ, dofCorrect = lin.N / lin.dof)
end


#=Simple test function to verify the QR-related optimizations
Output should be identical to the standard methods
IN: Independent variable X matrix
OUT: Returns the covariance matrix=#
getHomoskedΣSlow(X::Matrix{Float64}, ε::Vector{Float64})::Matrix{Float64} =
  (X'*X)\Matrix{Float64}(I,size(X,2),size(X,2))*(ε⋅ε)/(size(X,1)-size(X,2))

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
getHomoskedΣSlow(lin::FMLM)::Matrix{Float64} =
  getHomoskedΣSlow(lin.X, lin.Y.-lin.X * lin.β)


#########################getNeweyWest ##########OLS Only
#NOTE: To match Sandwich NeweyWest in R, set prewhite = FALSE, adjust = TRUE
function getNeweyWest!(X::Matrix{Float64}, xqr::FMQR, ε::Vector{Float64}, lag::Int,
  Σ::Matrix{Float64} = Matrix{Float64}(undef, xqr.K, xqr.K), dofCorrect::Float64 = 1.0)::Matrix{Float64}

  #pre-allocate for the spectral matrix

  Rv::Matrix{Float64} = Matrix{Float64}(undef, xqr.K, xqr.K) #pre-allocate working matrix
  RRInv::Matrix{Float64} = BLAS.gemm('N', 'T', xqr.RInv, xqr.RInv) #this is equivelent to [X'X]^-1

  #need to multiply through by the error
  Xe::Matrix{Float64} = X .* ε
  ST::Matrix{Float64} = BLAS.gemm('T','N',1.0/xqr.N, Xe, Xe)
  for v::Int ∈ 1:lag
    #overwrites Rv with (1/N)R'R
    BLAS.gemm!('T', 'N', 1.0/xqr.N, view(Xe, (v+1):(xqr.N), :),view(Xe, 1:(xqr.N-v), :), 0.0, Rv)
    #Rv .= view(Xe, (v+1):(xqr.N), :)' * view(Xe, 1:(xqr.N-v), :) .* (1.0/xqr.N)
    ST .+= (lag + 1 .- v)/(lag+1.) .* (Rv .+ Rv')
  end

  #this is [X'X]^-1S=[R'R]^-1S
  RRInvS::Matrix{Float64} = BLAS.gemm('N', 'N', RRInv, ST)

  #finally we have T[X'X]^-1S[X'X]^-1
  BLAS.gemm!('N','N',Float64(xqr.N), RRInvS, RRInv, 0.0, Σ)
  #println("$(diag(Σ .* dofCorrect))")
  return Σ .* dofCorrect
end


function getNeweyWest!(lin::FMLM, lag::Int,
  Σ::Matrix{Float64} = Matrix{Float64}(undef, lin.K, lin.K))::Matrix{Float64}

  return getNeweyWest!(lin.X, lin.xqr, lin.ε, lag, Σ, lin.N/lin.dof)
end

#helper function in case the input requires a single argument function
getNeweyWestFunc(lag::Int) = (lin::FMLM)-> getNeweyWest!(lin, lag)

#testing for NeweyWest
#NOTE: assumes dofcorrect value
function getNeweyWestSlow(X::Matrix{Float64}, ε::Vector{Float64}, lag::Int, dofCorrect::Float64)
  N::Int, K::Int = size(X)

  #first form the spectral matrix

  #initial value for spectral matrix
  xxt::Matrix{Float64} = Matrix{Float64}(undef, K,K)
  xxt .= zeros(K,K)
  for t::Int ∈ 1:N
    xxt .+= (X[t, :] * X[t, :]') .* ε[t]^2
  end
  ST::Matrix{Float64} =  xxt ./ N

  #make the rest of the spectral matrix
  for v::Int ∈ 1:lag
    xxt .= zeros(K,K)
    for t::Int ∈ (1+v):(N)
      xxt .+= (X[t, :] * X[t-v, :]') .* ε[t]*ε[t-v]
    end
    xxt ./= N
    ST .+= (lag + 1. - v)/(lag + 1.) .* (xxt .+ xxt')
  end

  XXInv::Matrix{Float64} = (X' * X) \ Matrix{Float64}(I,K,K)

  ret::Matrix = N*XXInv*ST* XXInv
  return ret * dofCorrect
end

getNeweyWestSlow(lin::FMLM, lag::Int) = getNeweyWestSlow(lin.X, lin.Y .- lin.X * lin.β, lag, lin.N/lin.dof)

getNeweyWestFuncSlow(lag::Int)::Function = (lin::FMLM)-> getNeweyWestSlow(lin, lag)
###################getWhiteERRORs ###OLS only

function getWhiteΣ!(xqr::FMQR, ε::Vector{Float64},
  Σ::Matrix{Float64} = Matrix{Float64}(undef, xqr.K, xqr.K), dofCorrect::Float64=1.0)::Matrix{Float64}

  Λ::Vector{Float64} = ε.^2.0

  QRtInv::Matrix{Float64} = BLAS.gemm('N','T',xqr.Q,xqr.RInv) #Q(R')^-1
  RInvQΛ::Matrix{Float64} = QRtInv'

  #loop to scale the matrix by Λ (modified part of modified white)
  @fastmath for j ∈ 1:xqr.N, i∈1:xqr.K
    RInvQΛ[i,j] *= Λ[j]
  end

  #final multiplicaiton and assignment
  return BLAS.gemm!('N','N',1.0,RInvQΛ, QRtInv,0.0,Σ) .* dofCorrect #R^-1Q'ΛQ(R')^-1
end

function getWhiteΣ!(lin::FMLM,
  Σ::Matrix{Float64} = Matrix{Float64}(undef, lin.K, lin.K))::Matrix{Float64}

  return getWhiteΣ!(lin.xqr, lin.ε, Σ, lin.N/lin.dof)
end

#for testing purposes only
function getWhiteΣSlow(X::Matrix{Float64}, ε::Vector{Float64})
  XXInv::Matrix{Float64} = (X' * X)\Matrix{Float64}(I, size(X,2), size(X,2))

  return XXInv * X' * diagm(ε) .^ 2.0 * X * XXInv
end

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
getWhiteΣSlow(lin::FMLM) = getWhiteΣSlow(lin.X, lin.Y.-lin.X * lin.β)

###########################getClustered SEs
#method for getting clustered errors
#from http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf
#forumala: [X'X]^-1*B*[X'X]^-1 where B= sum over G (Xg'*εg*εg'*Xg) and g is indexed for G clusters
function getClustered!(X::Matrix{Float64}, xqr::FMQR, ε::Vector{Float64}, clusters::Vector{C},
    Σ::Matrix{Float64} = Matrix{Float64}(undef, xqr.K, xqr.K),
    dof::Int = xqr.N-xqr.K)::Matrix{Float64} where C<:FMData

  #convenience vectors and values
  cVars::Vector{C} = unique(clusters)
  G::Int = length(cVars)
  cTable::Dict = Dict(cVars[i] => i for i::Int ∈ 1:G)
  clusterCode::Vector{Int} = ((x::C)->cTable[x]).(clusters)

  #iterate through the groups
  B::Matrix{Float64} = zeros(xqr.K, xqr.K)
  @fastmath for i::Int ∈ 1:G
    Xg::SubArray{Float64,2} = view(X, clusterCode .== i, :)
    ug::SubArray{Float64,1} = view(ε, clusterCode .== i)
    B .+= Xg'*ug*ug'*Xg
  end

  # Calc [X'X]^-1*B*[X'X]^-1
  RRInv::Matrix{Float64} =  BLAS.gemm('N','T',xqr.RInv,xqr.RInv)
  RRInvB::Matrix{Float64} = BLAS.gemm('N','N',RRInv,B)

  dofCorrect::Float64 = G/(G-1.)*(xqr.N-1.)/dof #small sample correction
  BLAS.gemm!('N','N',dofCorrect,RRInvB,RRInv,0.0,Σ)

  @fastmath for i ∈ 1:xqr.K
    Σ[i,i] = abs(Σ[i,i]) > 10. ^ -20. ? Σ[i,i] : 0.0
    if Σ[i,i] < 0.
      error("Negative variance ($i): $(diag(Σ))")
    end
  end

  #println("$(diag(Σ))")
  return Σ
end

function getClustered!(lin::FMLM, Σ::Matrix{Float64} = Matrix{Float64}(undef, lin.K, lin.K))::Matrix{Float64}

  return getClustered!(lin.X, lin.xqr, lin.ε, lin.clusters, Σ, lin.dof)
end

###################getMWErrors!##############################
#Main method to get the modified white SEs Does #[X'X]^-1X'ΛX[X'X]^-1


#=Main method to calculate the modified white SEs
IN: A FMQR decomposition object, a residuals vector,
optional memory for a covariance matrix
OUT: Writes and returns the covariance matrix=#
function getModWhiteΣ!(xqr::FMQR, ε::Vector{Float64},
  Σ::Matrix{Float64} = Matrix{Float64}(undef, xqr.K, xqr.K), dofCorrect::Float64 = 1.0)::Matrix{Float64}

  Λ::Vector{Float64} = similar(ε)

  project!(xqr,Λ) #get the projection diagonal
  Λ .= (ε./(1.0.-Λ)).^2.0

  QRtInv::Matrix{Float64} = BLAS.gemm('N','T',xqr.Q,xqr.RInv) #Q(R')^-1
  RInvQΛ::Matrix{Float64} = QRtInv'

  #loop to scale the matrix by Λ (modified part of modified white)
  @fastmath for j ∈ 1:xqr.N, i∈1:xqr.K
    RInvQΛ[i,j] *= Λ[j]
  end

  #final multiplicaiton and assignment
  return BLAS.gemm!('N','N',dofCorrect,RInvQΛ, QRtInv,0.0,Σ) #R^-1Q'ΛQ(R')^-1
end

#=helper method which extracts the required components from a linear model
IN: A linear model, optionally memory for the covariance matrix
OUT: Writes and allocates the covariance matrix=#
function getModWhiteΣ!(lin::FMLM,
  Σ::Matrix{Float64} = Matrix{Float64}(undef, lin.K, lin.K))::Matrix{Float64}

  return getModWhiteΣ!(lin.xqr, lin.ε, Σ, lin.N/lin.dof)
end


#=Simple test function to verify the QR-related optimizations
Output should be identical to the standard methods
IN: Independent variable X matrix
OUT: Returns the covariance matrix=#
function getModWhiteΣSlow(X::Matrix{Float64}, ε::Vector{Float64})
  XXInv::Matrix{Float64} = (X' * X)\Matrix{Float64}(I,size(X,2),size(X,2))

  P::Matrix{Float64} = X * (XXInv) * X'
  return XXInv * X' * diagm(ε  ./ (1.0 .- diag(P))) .^ 2.0 * X * XXInv

end

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
getModWhiteΣSlow(lin::FMLM) = getModWhiteΣSlow(lin.X, lin.Y.-lin.X * lin.β)


  #################testing methods for this module
  #NOTE: These need to be re-written
function LMtest()

  #Number of Observations
  N::Int = 200
  K::Int = 2

  #Allocate
  X::Matrix{Float64} = Matrix{Float64}(undef, N,K)
  Y::Vector{Float64} = Vector{Float64}(undef, N)
  ε::Vector{Float64} = similar(Y)

  #parameters for the simulation
  σ2::Vector{Float64} = [2.0^2.0 for i::Int ∈ 1:N]
  σ2[1:ceil(Int, N/2)] .= 0.5^2.0
  β::Vector{Float64} = [(1.0 * i) for i::Int ∈ 1:K]
  #println("β: ", β)
  #println("σ2: $σ2")

  #this will hold the sampled beta
  e::Vector{Float64} = similar(ε)

  #run the simulation
  X[:,1] .= 1.0 #intercept
  for i ∈ 2:K
    X[:,i] .= rand(Uniform(),N)
  end

  ε .= map((s2::Float64)->rand(Normal(0.0,s2^0.5)),σ2)
  Y .=  X*β .+ ε

  #get the linear model
  lin::FMLM = FMLM(X, Y)

  #get the homoskedastic SEs
  ΣHomosked::Matrix{Float64} = getHomoskedΣ!(lin)
  ΣHomoskedSlow::Matrix{Float64} = getHomoskedΣSlow(lin)

  #print the coefficients
  println("Coefficients: ",lin.β)
  println("Homoskedastic Errors: ", diag(ΣHomosked).^.5)
  println("Check: ", diag(ΣHomoskedSlow).^.5)

  #get the modified white SEs
  ΣWhite::Matrix{Float64} = similar(ΣHomosked)
  getModWhiteΣ!(lin, ΣWhite)
  ΣWhiteSlow::Matrix{Float64} = getModWhiteΣSlow(lin)

  #print the coefficients
  println("Modified White Errors: ",diag(ΣWhite).^.5)
  println("Check: ", diag(ΣWhiteSlow).^.5)

    #=Π1::Matrix{Float64} = Matrix{Float64}(K,2)
  getCoeff!(lin.xqr, [Y 2.*Y], Π1)
  println("1st stage Z: $(lin.X)")
  println("1st Stage X: $([Y 2.*Y])" )
  println("1st Stage Test Coef: $Π1" )
  Ξ::Matrix{Float64} = getResid!(lin.xqr, [Y 2.*Y], Π1)
  println("1st Stage Resid: $Ξ" )=#


  #test the project routines
  Pa::Vector{Float64} = Vector{Float64}(undef, N)
  PM::Matrix{Float64} = Matrix{Float64}(undef, N,N)
  PS::Matrix{Float64} = Matrix{Float64}(undef, N,N)

  project!(X,PM)
  project!(X,Pa)

  projectSlow!(X,PS)
  println("P: ", Pa[1:5])
  println("P (from full matrix): ", PM[1:3,1:3])
  println("P Slow: ", diag(PS)[1:10])
end



if DEBUG_FMMOD
  #LMtest()
  IVTest()
  #IOTest()
  #print(perfTestStr(10^6))
end
