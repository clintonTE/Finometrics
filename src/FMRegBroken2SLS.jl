#functions and extensiosn for 2SLS
#WARNING: Probably broken, hasn't been maintained in a while

######2SLS constructor
#Holds the modelling componenets for a 2SLS model (except errors)
#Variable descriptions are in the struct
# Runs standard 2SLS using QR optimization
mutable struct FM2SLS <: FMModel
    zaqr::FMQR #QR decomposition of Za (Z|W)
    xaqr::FMQR #QR decomposition of fitted 1st stage values (X̃|W)

    X::Matrix{Float64} #Endogeneous covariates
    W::Matrix{Float64} #Exogeneous covariates
    Y::Vector{Float64} #Dependent var
    Z::Matrix{Float64} #IV data

    N::Int #number of data points
    K::Int # of endogeneous covariates
    L::Int # of instrumental variables
    KW::Int # of exogeneous covariates
    dof::Int

    Π1::Matrix{Float64} #first stage results
    Ξ1::Matrix{Float64} #first stage residuals

    δ2::Vector{Float64} #2SLS coefficients
    ξ2::Vector{Float64} #2SLS residuals

    XNames::Vector{Symbol} #names of the X variables
    WNames::Vector{Symbol} #names of the W variables
    YName::Symbol #names of the Y variables
    ZNames::Vector{Symbol} #names of the Z variables

    #FM2SLS()
end


#= Convenience constructor for FM2SLS that allocates a 0 element array for the
case of no exogeneous covariates. NOTE: no exog covariates implies no intercept,
so there probably should be  some exogeneous covariates
#IN: X Matrix, Y Vector, and the Z matrix of Instrumental variables
#OUT: FM2SLS object=#
FM2SLS(X::Matrix{Float64},Y::Vector{Float64},Z::Matrix{Float64};
    XNames::Vector{Symbol} = Symbol.(:X, 1:size(X,2)),
    YName::Symbol = :Y,
    ZNames::Vector{Symbol} = Symbol.(:Z, 1:size(Z,2)))::FM2SLS =
        FM2SLS(X,Matrix{Float64}(undef, Size(X,1),0),Y,Z)

#=Main constructor for FM2SLS object
IN: Endogeneous matrix X, Exogeneous matrix W, Dependent vector Y,
instrumental zariable matrix Z
OUT: A 2SLS object

#NOTE: not sure about the dof
=#
function FM2SLS(X::Matrix{Float64}, W::Matrix{Float64},
  Y::Vector{Float64},Z::Matrix{Float64};
  XNames::Vector{Symbol} = Symbol.(:X, 1:size(X,2)),
  WNames::Vector{Symbol} = Symbol.(:W, 1:size(W,2)),
  YName::Symbol = :Y,
  ZNames::Vector{Symbol} = Symbol.(:Z, 1:size(Z,2)),
  dof::Int = size(X,1) - size(X,2) - size(W,2))::FM2SLS

  #memory pre-allocaiton
  N::Int = size(X,1) #number of data points
  K::Int = size(X,2) # of endogeneous variables
  L::Int = size(Z,2) # of instrumental variables
  KW::Int = size(W,2)# of exogeneous covariates

  Π1::Matrix{Float64} = Matrix{Float64}(undef, L+KW,K)
  X̂::Matrix{Float64} = Matrix{Float64}(N,K) #holds fitted first stage values
  Ξ1::Matrix{Float64} = similar(X̂) #holds first stage residuals
  δ2::Vector{Float64} = Vector{Float64}(undef, K+KW)
  ξ2::Vector{Float64} = similar(Y)
  #Xa::Matrix{Float64} = Matrix{Float64}(N,K+KW) #fitted 1st stage plus exog covariates

  #start with setting up and running the first stage
  ZW::Matrix{Float64} = [Z W]
  zaqr::FMQR = FMQR(ZW, keepX=false) #augment the matrices to ZW and get the QR
  getCoeff!(zaqr, X, Π1) #get the first stage coefficients
  BLAS.gemm!('N','N',1.0,ZW,Π1,0.0,X̂) #get the fitted first stage values
  getResid!(ZW, X, Π1, Ξ1) #get the residuals

  #now run the second stage
  xaqr::FMQR = FMQR([X̂ W], keepX = false) #augment X̂W and get the QR
  getCoeff!(xaqr, Y, δ2) #get the 2SLS coefficients
  getResid!([X W], Y, δ2, ξ2) #get the 2SLS residuals

  #Finally, return the constructed model
  return FM2SLS(zaqr, xaqr, X, W, Y, Z, N, K, L, KW, dof, Π1, Ξ1, δ2, ξ2,
    XNames, WNames, YName, ZNames)
end




#Constructor and helper method which gets required info from DataFrame
#IN: The source dataframe, a formula expression/symbol for X, W, Y, and Z
# the dependent variable as a symbol
#OUT: A FMLM object
function FM2SLS(df::DataFrame, XExpr::T, WExpr::U, YSym::Symbol, ZExpr::V;
    XNames::Vector{Symbol} = Vector{Symbol}(),
    WNames::Vector{Symbol} = Vector{Symbol}([:Intercept]),
    YName::Symbol = :Y,
    ZNames::Vector{Symbol} = Vector{Symbol}(),
    containsmissings=true)::FM2SLS where {T <: FMExpr, U <: FMExpr, V <: FMExpr}

    #get the list of symbols for X, W and Z
    XSym::Vector{Symbol} = getSymbolsFromExpr(XExpr)
    WSym::Vector{Symbol} = getSymbolsFromExpr(WExpr)
    ZSym::Vector{Symbol} = getSymbolsFromExpr(ZExpr)

    #make a view so if we drop nulls it doesn't affect the original
    dfSub::SubDataFrame = view(df[ [XSym; WSym; YSym; ZSym]],1:size(df,1),:)

    # if we are going to drop nulls
    if containsmissings
        dfSub = cloakmissings(dfSub)
    end

    #get the model matrices, which will expand all factors
    XModelMatrix::Matrix{Float64} = getModelMatrix(dfSub, XExpr)
    WModelMatrix::Matrix{Float64} = getModelMatrix(dfSub, WExpr)
    ZModelMatrix::Matrix{Float64} = getModelMatrix(dfSub, ZExpr)

    #this is a check to make sure the factor expansion follows other columns
    #We must repeat the check on each variable matrix
    for symVec::Vector{Symbol} ∈ [XSym, WSym, ZSym]
        if length(symVec) > 1
            for i ∈ 2:length(symVec) #all types
                if typeof(dfSub[symVec[i-1]])<:CategoricalArray &&
                    !(typeof(dfSub[symVec[i]])<:CategoricalArray)
                    @warn "Type of $(typeof(dfSub[i])) is preceeded by type of
                        factor. Column names are likely MISALIGNED. Put factors
                        after numerical variables in RHS expressions."
                end
            end
        end
    end

    #now create the namesfor the expanded matrix
    XNamesFull::Vector{Symbol} = Vector{Symbol}(undef, size(XModelMatrix,2))
    for i ∈ 1:length(XNamesFull)
        XNamesFull[i] = i≤length(XNames) ? XNames[i] : Symbol(:X,i)
    end

    WNamesFull::Vector{Symbol} = Vector{Symbol}(undef, size(WModelMatrix,2))
    for i ∈ 1:length(WNamesFull)
        WNamesFull[i] = i≤length(WNames) ? WNames[i] : Symbol(:W,i)
    end

    ZNamesFull::Vector{Symbol} = Vector{Symbol}(undef, size(ZModelMatrix,2))
    for i ∈ 1:length(ZNamesFull)
        ZNamesFull[i] = i ≤ length(ZNames) ? ZNames[i] : Symbol(:Z,i)
    end

    #get the appropriate matrices and construct the FM2SLS object
    return FM2SLS(XModelMatrix, WModelMatrix, Vector{Float64}(undef, dfSub[YSym]),
        ZModelMatrix, XNames=XNamesFull, WNames=WNamesFull,
        YName=YName, ZNames=ZNamesFull)

end

#default constructor to help with memory management
FM2SLS() = FM2SLS(FMQR(), FMQR(), Matrix{Float64}(), Matrix{Float64}(),
    Vector{Float64}(), Matrix{Float64}(), 0, 0, 0, 0, Matrix{Float64}(),
    Matrix{Float64}(), Vector{Float64}(), Vector{Float64}(),
    Vector{Symbol}(), Vector{Symbol}(), Symbol(), Vector{Symbol}())::FM2SLS

###########################2SLS Helper functions

#=function get1stStage
#The purpose is to get the first stage of the IV as a linear model
#Unfortunately, we need to copy some of the data, so this can be expensive,
#although Z and W should only be copied once
#IN: an IV object
#OUT: a vector of FMLM objects, one for each focal (X) variable in the
#initial regression=#
function get1stStage(iv::FM2SLS)::Vector{FMLM}
    ZW::Matrix{Float64} =  [iv.Z iv.W]
    iv1st::Vector = Vector{FMLM}()

    for i ∈ 1:length(iv.XNames)
        push!(iv1st, FMLM(iv.zaqr, [iv.Z iv.W], iv.X[:,i],
            iv.zaqr.N, iv.zaqr.K, iv.Π1[:,i], iv.Ξ1[:,i],
            [iv.ZNames; iv.WNames], iv.XNames[i]))
    end

    return iv1st
end

#=gets the residuals in the 2SLS framework
Note this is a wrapper function which can accept multiple independent
and multiple dependent variables. NOTE the change in the naming convention
relative to other verisons of this method (due to the 2SLS context)
IN: QR decomposition of the RHS variables, LHS dependent variable matrix,
a coefficient matrix, optionally the memory for the residual matrix
OUT: Writes and returns the residual matrix=#
function getResid!(ZW::Matrix{Float64}, X::Matrix{Float64}, Π1::Matrix{Float64},
  Ξ1::Matrix{Float64}=similar(X))::Matrix{Float64}

  #get the residuals by reference
  for i ∈ 1:size(X,2) #get the vector of reisduals multiple times
    getResid!(ZW, view(X,:,i), view(Π1,:,i), view(Ξ1,:,i))
  end
  return Ξ1
end



#=This is a wrapper function which gets a matrix of coefficients for use in
# the first stage of 2SLS. NOTE the change in the naming convention
relative to other verisons of this method (due to the 2SLS context)
IN: QR decomposition of the RHS variables, LHS dependent variable matrix,
optionally the memory for the coefficient matrix
OUT: Writes and returns the coefficient matrix=#
function getCoeff!(zaqr::FMQR, X::Matrix{Float64},
  Π1::Matrix{Float64} = Matrix{Float64}(undef, zaqr.K, size(X,2)))::Matrix{Float64}

  for i ∈ 1:size(X,2) # for each dependent variable vector
    getCoeff!(zaqr, view(X,:,i), view(Π1,:,i)) #run the regression as needed
  end

  return Π1
end

################################2SLS Standard Errors

#Convenience method for above: IN: A FM2SLS object #OUT: A covariance matrix
getHomoskedΣSlow(iv::FM2SLS)::Matrix{Float64} =
  getHomoskedΣSlow([[iv.Z iv.W]*iv.Π1 iv.W], iv.Y.-[iv.X iv.W] * iv.δ2)

#=helper method for the above function which extracts the QR decomposition
from a FM2SLS object
IN: A linear model, optionally memory for the covariance matrix
OUT: WRites and returns the covariance matrix=#
function getHomoskedΣ!(iv::FM2SLS,
  Σ::Matrix{Float64}=Matrix{Float64}(undef, iv.xaqr.K, iv.xaqr.K))::Matrix{Float64}

  return getHomoskedΣ!(iv.xaqr, iv.ξ2, Σ, iv.N / iv.dof)
end


#=helper method which extracts the required components from a 2SLS model
IN: A 2SLS model, optionally memory for the covariance matrix
OUT: Writes and allocates the covariance matrix=#
function getModWhiteΣ!(iv::FM2SLS,
  Σ::Matrix{Float64}=Matrix{Float64}(undef, iv.xaqr.K, iv.xaqr.K))::Matrix{Float64}

  return getModWhiteΣ!(iv.xaqr, iv.ξ2, Σ, iv.N/iv.dof)
end


#Convenience method for above: IN: A 2SLS object #OUT: A covariance matrix
getModWhiteΣSlow(iv::FM2SLS)::Matrix{Float64} =
  getModWhiteΣSlow([[iv.Z iv.W]*iv.Π1 iv.W], iv.Y.-[iv.X iv.W] * iv.δ2)


###########################IO
#if we get serious about 2SLS again these can go into the IO file

function texTable(models::Vector{FM2SLS}, getΣ::T where T<:Union{Function, Vector{Function}},
    rows::Vector{Symbol};
    titleCaption::String = "",
    caption::String = "",
    colNames::Vector{Vector{String}} = [["C$i" for i ∈ 1:length(models)]],
    contentRowNames::Vector{String} = String.(rows),
    descRowNames::Vector{String}=Vector{String}(),
    descContent::Vector{Vector{String}}=Vector{Vector{String}}(),
    notes::Vector{String} = Vector{String}(),
    arrayStretch::Float64 = 1.5,
    lineSpacer::String = "\\\\",
    summaryMathMode::Bool = true,
    widthColNames::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(colNames[i])),1:length(colNames)),
    alignmentColNames::Vector{Vector{String}} = #contains the number of columns for each entry
      broadcast((i::Int)->["r" for i ∈ 1:length(colNames[i])],1:length(colNames)),
    widthDescContent::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(descContent[i])),1:length(descContent)),
    stars::Bool=true,
    starLvls::Vector{Float64} = [.9, .95, .99],
    starLegend::String = stars ? DEFAULT_STAR_LEGEND : "",
    starStrings::Vector{String} =
      ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    scaling::Vector{Float64}=ones(length(rows)),
    decimalDigits::Int = 2,
    colHeaderName::Vector{String} = ["" for i ∈ 1:length(colNames)],
    nakedTable=nakedTable,
    alignmentstring::String = string(" l | ", join(["r" for i ∈ 1:length(descContent[1])])))

    numCols = length(models)
    numContentRows = length(rows)

    #Pre-allocate the vector of SE errors
    modelsσ::Vector{Vector{Float64}} =
    [Vector{Float64}(models[i].K+models[i].KW) for i∈1:numCols]


    #pull out the β coefficients and N
    modelsβ::Vector{Vector{Float64}} = [models[i].δ2 for i ∈ 1:numCols]
    modelsN::Vector{Int} = [models[i].N for i ∈ 1:numCols]
    modelsXWNames::Vector{Vector{Symbol}} = [[models[i].XNames; models[i].WNames]  for i ∈ 1:numCols]

    #get the standard errors
    for c ∈ 1:numCols
        if typeof(getΣ) <: Function
          modelsσ[c] .= sqrt.(diag(getΣ(models[c])))
        else
          modelsσ[c] .= sqrt.(diag(getΣ[c](models[c])))
        end
    end

    if length(starLegend) > 0
        notes = [starLegend; notes]
    end

    #get the content matrices
    content::Vector{Matrix{String}} = getContentMatrices!( modelsβ, modelsσ, modelsXWNames, modelsN, rows,
      stars=stars, starLvls=starLvls, scaling=scaling,
      decimalDigits=decimalDigits, starStrings=starStrings)

    texTable(titleCaption = titleCaption,caption=caption, colNames=colNames,
      contentRowNames=contentRowNames, content=content, descRowNames=descRowNames,
      descContent=descContent, notes=notes, arrayStretch=arrayStretch, lineSpacer=lineSpacer,
      summaryMathMode=summaryMathMode, widthColNames=widthColNames,
      widthDescContent=widthDescContent, colHeaderName=colHeaderName,
      alignmentColNames=alignmentColNames, nakedTable=nakedTable, alignmentstring=alignmentstring)
end


function IVTest()

  srand(11)

  N::Int = 200 #Number of Observations
  K::Int = 3 # of endogneous covariates
  L::Int = 4 # of instruments
  KW::Int = 5# of exog covariates

  #Allocate
  X::Matrix{Float64} = Matrix{Float64}(undef, N,K)
  Z::Matrix{Float64} = Matrix{Float64}(undef, N,L)
  W::Matrix{Float64} = Matrix{Float64}(undef, N,KW)
  Y::Vector{Float64} = Vector{Float64}(undef, N)
  A::Vector{Float64} = Vector{Float64}(undef, N)

  Π1::Matrix{Float64} = Matrix{Float64}(undef, L+KW,K)
  π2::Vector{Float64} = Vector{Float64}(undef, L+KW)
  Ξ1::Matrix{Float64} = similar(X)
  e2::Vector{Float64} = similar(Y)

  #parameters for the simulation
  σ2::Vector{Float64} = [2.0^2 for i::Int ∈ 1:N]
  σ2[1:ceil(Int, N/2)] .= 0.5^2.0

  #generate π2 coefficients (pattern is 1 times the index, ie 1 2 3 ...)
  π2 .= [(10*i) for i::Int ∈ 1:(L+KW)]
  π2[L+1] .= 1.0 #intercept term
  #println("π2: ",π2)

  #generate Π1, pattern is δ2 times 1/10 the index, ie ie .1δ2 .2δ2 .3δ2 ...)
  for i ∈ 1:(L+KW), j∈1:K
    Π1[i,j] .= i + j^2
  end
  Π1[L+1,:] .= 1.0 #set the intercept term
  #display(Π1)

  #generate the exogeneous variables
  for i ∈ 2:KW
    W[:,i] .= rand(Uniform(-10,10),N)
  end
  W[:,1] .= 1.0 #intercept

  #generate the instrumental variables
  for i ∈ 1:L
    Z[:,i] .= rand(Uniform(-10,10),N)
  end

  #map((s2::Float64)->rand(Normal(0.0,s2^0.5)),σ2)
  #generate the errors
  for i ∈ 1:N
    for j ∈ 1:K
      Ξ1[i,j] = rand(Normal(0.0,σ2[i]^0.5))
    end
    e2[i] = rand(Normal(0.0,σ2[i]^0.5))
  end

  #println(Ξ1)
  #populate Y and X
  X .= [Z W]*Π1 .+ Ξ1
  Y .= [Z W]*π2 .+ e2

  #add some endogeneity
  for i ∈ 1:K
    A .= rand(Uniform(-i,i),N)
    X[:,i] .+= A
    Y .+= A .* 10
  end


  #probably is a more graceful way to do this
  XNames::Vector{Symbol} = [Symbol("X$i") for i = 1:K]
  WNames::Vector{Symbol} = [Symbol("W$i") for i = 1:KW]
  ZNames::Vector{Symbol} = [Symbol("Z$i") for i = 1:L]


  iv::FM2SLS = FM2SLS(X,W,Y,Z, XNames=XNames, YName=:Y, WNames=WNames, ZNames=ZNames)
  println("coefficients: \n", [XNames; WNames] ,"\n", iv.δ2)

  #get the homoskedastic SEs
  ΣHomosked::Matrix{Float64} = getHomoskedΣ!(iv)
  ΣHomoskedSlow::Matrix{Float64} = getHomoskedΣSlow(iv)

  #print the coefficients
  println("Homoskedastic Errors: ", diag(ΣHomosked).^.5)
  println("Check: ", diag(ΣHomoskedSlow).^.5)

  #get the modified white SEs
  ΣWhite::Matrix{Float64} = getModWhiteΣ!(iv)
  ΣWhiteSlow::Matrix{Float64} = getModWhiteΣSlow(iv)

  #print the coefficients
  println("Modified White Errors: ",diag(ΣWhite).^.5)
  println("Check: ", diag(ΣWhiteSlow).^.5)

  iv1st::Vector{FMLM} = get1stStage(iv)
  println("1st Stage results: \n $(iv1st[1].XNames)")
  for iv1 ∈ iv1st
      println("coefficients $(iv1.YName):  $(iv1.β )")
      println("Error (Homosked $(iv1.YName)): $(diag(getHomoskedΣ!(iv1)).^.5)")
  end



  oStream::IOStream = open(
    "C:\\Users\\Clinton\\Dropbox\\Projects\\Summer17RTester\\Summer17RTester\\IVTestOut.csv","w+")


  write(oStream, [string.(iv.XNames); string.(iv.WNames)]::Vector{String})
  writecsv(oStream, [X W Y Z])
  close(oStream)
end
