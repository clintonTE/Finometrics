using Revise, Dates, DataFrames, GLM
#NSymbol = Union{Nothing, Symbol}
#NInt = Union{Nothing, Int} #NOTE: consider deleting this later
import Base: +, -, ==, >, <, ≥, ≤, length, isless, isequal
import StatsModels: implicit_intercept


#include("Finometrics.jl")
using Revise
const FMExpr = Union{Symbol,Expr,Nothing}

macro mpar(cond, expr)
  quote
    if $(esc(cond))
        :($(Threads.@threads($expr)))
    else
        :($($expr))
    end
  end
end

#  StatsModels.implicit_intercept(::Type{<:Any}) = true
#=function getModelMatrix(df::T, f::FormulaTerm)::Matrix{Float64} where
  T <: AbstractDataFrame

  mf::ModelFrame = ModelFrame(f, df)
  #println("frame created")
  return ModelMatrix(mf).m
end

#Same as above but allows for an expression
function getModelMatrix(df::T, rhs::V)::Matrix{Float64} where
  {T <: AbstractDataFrame, V <: FMExpr}

  #special case which crashes Formula
  if exp == Symbol("")
    return ones(Float64,size(df,1),1)
  end

  #WARNING: This is a HACK! Replace with `nothing` or something else (See StatsModels github)
  #or just redo the whole thing in the GLM framework
  lhs_hack::Symbol = names(df)[1]
  f::FormulaTerm = @eval(@formula($lhs_hack~ $rhs))

  return getModelMatrix(df, f)
end

function getModelMatrix(df::T, f::FormulaTerm)::Matrix{Float64} where
  T <: AbstractDataFrame
  f = apply_schema(f, schema(f,df), StatisticalModel)
  m = modelcols(f.rhs, df)
  return m
end

#Same as above but allows for an expression
function getModelMatrix(df::T, rhs::V)::Matrix{Float64} where
  {T <: AbstractDataFrame, V <: FMExpr}

  #special case which crashes Formula
  if rhs == Symbol("")
    return ones(Float64,size(df,1),1)
  end

  #WARNING: This is a HACK! Replace with `nothing` or something else (See StatsModels github)
  #or just redo the whole thing in the GLM framework
  #Will be able to replace the LHS w/0 in subsequent versions
  lhs_hack::Symbol = names(df)[1]
  f::FormulaTerm = @eval(@formula($lhs_hack ~ $rhs))

  return getModelMatrix(df, f)
end


function testMM(df::DataFrame, loops=1)
  #println(df)
  x = getModelMatrix(df, :(x+y+f))

  (loops==1) && display(x)

  return sum(x)
end

function testMM(loops::Int; N::Int=1000, G::Int=50)
  local tot::Float64 = 0.0

  df = DataFrame(x = 1001:(1000+N), y = rand(N), z=rand(N),
    f = (i->Symbol(:q,i÷(N ÷ G))).(1:N))

  @time for i ∈ 1:loops
    tot+=testMM(df, loops)/N/G
    df.y .+= rand()
  end

  println(tot)
end

testMM(10_000, N=1000, G=10)=#




#=using Finometrics
using Revise


#use this to test the lag functions
function testlaganddifference(N = 10_000)
  rowindex = collect(1:N)
  periodcol = deepcopy(rowindex)


  #make a test frame
  df = DataFrame(rowindex=rowindex, period = periodcol, group1 = periodcol .÷ 17,
    group2 = periodcol .÷ 29, target1 = rand(N), target2 = rand(N))

  #sort!(df, [:group1, :period, :group2])
  #CASE 1: create lagged field Lgroup1
  lagwithin!(df, :target1, :group1, :period)

  #CASE 2: create lagged fields L2target1 and L2target2
  lagwithin!(df, [:target1, :target2], [:group1, :group2], :period, lags=2, sorted=false)

  #CASE 3: difference :Ltarget1 and create DLtarget1 and LLtarget1
  differencewithin!(df, :Ltarget1, :group1, :period, deletelag=false, sorted=true,
    laggedname=:LLtarget1e, differencedname=:DLtarget1e)

  #CASE 4: difference [:L2target1, :L2target2] creating
  # [:DL2target1e, :DL2target2e], [:LL2target1e, :LL2target2e]
  differencewithin!(df, [:L2target1, :L2target2], [:group1, :group2], :period,
    deletelag=false, sorted=true, differencednames=[:DL2target1e, :DL2target2e],
    laggednames = [:LL2target1e, :LL2target2e])

  #validation for cases 1 and 3
  sort!(df, [:group1, :period])
  for subdf ∈ groupby(df, :group1)
    sN::Int = size(subdf, 1)
    for r ∈ 1:sN
      if r > 1 #CASE 1check
        @assert subdf[r-1, :target1] == subdf[r, :Ltarget1]
      else
        @assert ismissing(subdf[r, :Ltarget1])
      end

      if r > 2 #CASE 3 check
        @assert subdf[r-1, :Ltarget1] == subdf[r, :LLtarget1e]
        @assert subdf[r, :Ltarget1] - subdf[r, :LLtarget1e] == subdf[r, :DLtarget1e]
      else
        @assert ismissing(subdf[r, :LLtarget1e])
        @assert ismissing(subdf[r, :DLtarget1e])

      end
    end
  end

  #validation for cases 2 and 4
  sort!(df, [:group1, :group2, :period])
  for subdf ∈ groupby(df, [:group1, :group2])
    sN::Int = size(subdf, 1)
    for r ∈ 1:sN
      if r > 2 #CASE 2 check
        @assert subdf[r-2, :target1] == subdf[r, :L2target1]
        @assert subdf[r-2, :target2] == subdf[r, :L2target2]
      else
        @assert ismissing(subdf[r, :L2target1])
        @assert ismissing(subdf[r, :L2target2])
      end

      if r > 3 #CASE 4 check
        @assert subdf[r-1, :L2target1] == subdf[r, :LL2target1e]
        @assert subdf[r-1, :L2target2] == subdf[r, :LL2target2e]
        @assert subdf[r, :L2target1] - subdf[r, :LL2target1e] == subdf[r, :DL2target1e]
        @assert subdf[r, :L2target2] - subdf[r, :LL2target2e] == subdf[r, :DL2target2e]
      else
        @assert ismissing(subdf[r, :LL2target1e])
        @assert ismissing(subdf[r, :LL2target2e])
        @assert ismissing(subdf[r, :DL2target1e])
        @assert ismissing(subdf[r, :DL2target2e])

      end
    end
  end

  display(df)
end

#testlaganddifference()

function testYearQuarter()

  #test the first constructor
  yq1::YearQuarter = YearQuarter(1998,2)
  try
    YearQuarter(1985,7)
    @assert false
  catch err
    (typeof(err)<:AssertionError) && error("Error catching invalid values in first constructor")
  end

  #test the second constructor
  yq2::YearQuarter = YearQuarter(2000.4)
  try
    YearQuarter(2000.05)
    @assert false
  catch err
    (typeof(err)<:AssertionError) && error("Error catching invalid values in 2nd constructor")
  end

  #test addition of quarters
  for i ∈ 0:10
    @assert (yq1 + i) == (yq2 - (10-i))
  end
  #@assert yq1 + 10 == yq2
  #@assert yq2 - 10 == yq1

  @assert yq2 - yq1 == 10
  @assert yq1 < yq2
  @assert yq2 > yq1
  @assert yq2 ≥ yq1
  @assert yq1 ≥ yq1

  #println(yq1)
  #println(Float64(yq1))
end=#

function testYearMonth()

  #test the first constructor
  ym1::YearMonth = YearMonth(1998,2)
  try
    YearMonth(1985,14)
    @assert false
  catch err
    (typeof(err)<:AssertionError) && error("Error catching invalid values in first constructor")
  end

  #test the second constructor
  ym2::YearMonth = YearMonth(2008_04)
  try
    YearMonth(2000_14)
    @assert false
  catch err
    (typeof(err)<:AssertionError) && error("Error catching invalid values in 2nd constructor")
  end

  #test addition of quarters
  for i ∈ 0:122
    @assert (ym1 + i) == (ym2 - (122-i))
  end

  @assert ym2 - ym1 == 122
  @assert ym1 < ym2
  @assert ym2 > ym1
  @assert ym2 ≥ ym1
  @assert ym1 ≥ ym1

end

#@time for i ∈ 1:200_000
#  testYearMonth()
#end

function testymperformance(N::Int)
  ms::Vector{Int} = rand(1:12, N)
  ys::Vector{Int} = rand(0:20, N) .+ 1990
  addedmonths::Vector{Int} = rand(-120:120, N)
  local tot::Int
  local yms::Vector{YearMonth}
  local dts::Vector{Date}

  local addedMonths::Vector{Month} = (Month).(addedmonths)

  benchym::YearMonth = YearMonth(2000,1)
  benchdt::Date = Date(2000,1,1)

  @time begin
    yms = ((y::Int,m::Int)->YearMonth(y,m)).(ys,ms)
    yms .= yms .+ addedmonths
    tot = sum(yms .> benchym)
    print("\nYearMonth: $tot values found in")
  end

  @time begin
    dts = ((y::Int,m::Int)->Date(y,m,1)).(ys,ms)
    dts .= dts .+ addedMonths
    tot = sum(dts .> benchdt)
    print("Date: $tot values found in")
  end

end




#testymperformance(100_000_000)

include("Finometrics.jl")
using Distributions, LinearAlgebra, CuArrays

function LMtest(::Type{M}=Matrix{Float64}, ::Type{V}=Vector{Float64};
    N::Int = 200, K::Int = 2, testerrors::Bool = true,
    qrtype::Type = M) where {
    M<:AbstractMatrix, V<:AbstractVector}


  #Allocate
  X::M = M(undef, N,K)
  Y::V = V(undef, N)
  ε::V = similar(Y)

  #parameters for the simulation
  σ2::Vector{Float64} = [2.0^2.0 for i::Int ∈ 1:N]
  σ2[1:ceil(Int, N/2)] .= 0.5^2.0
  β::V = [(1.0 * i) for i::Int ∈ 1:K]
  #println("β: ", β)
  #println("σ2: $σ2")

  #this will hold the sampled beta
  e::Vector{Float64} = similar(ε)

  #run the simulation
  X[:,1] .= 1.0 #intercept
  for i ∈ 2:K
    X[:,i] .= V(rand(Uniform(),N))
  end

  ε = ((s2::Float64)->rand(Normal(0.0,s2^0.5))).(σ2)
  Y =  X*β .+ ε

  #get the linear model
  lin::Finometrics.FMLM = Finometrics.FMLM(X, Y, qrtype=qrtype)

  if testerrors
    #get the homoskedastic SEs
    ΣHomosked::Matrix{Float64} = Finometrics.getHomoskedΣ!(lin)
    ΣHomoskedSlow::Matrix{Float64} = Finometrics.getHomoskedΣSlow(lin)

    #print the coefficients
    println("Coefficients: ",lin.β)
    println("Homoskedastic Errors: ", diag(ΣHomosked).^.5)
    println("Check: ", diag(ΣHomoskedSlow).^.5)

    #get the modified white SEs
    ΣWhite::Matrix{Float64} = similar(ΣHomosked)
    Finometrics.getModWhiteΣ!(lin, ΣWhite)
    ΣWhiteSlow::Matrix{Float64} = Finometrics.getModWhiteΣSlow(lin)

    #print the coefficients
    println("Modified White Errors: ",diag(ΣWhite).^.5)
    println("Check: ", diag(ΣWhiteSlow).^.5)

    #get the nw SEs
    ΣNW::Matrix{Float64} = similar(ΣHomosked)
    Finometrics.getNeweyWest!(lin, 3, ΣNW)
    ΣNWSlow::Matrix{Float64} = Finometrics.getNeweyWestSlow(lin, 3)

    #print the coefficients
    println("NW Errors: ",diag(ΣNW).^.5)
    println("Check: ", diag(ΣNWSlow).^.5)



    #test the project routines
    Pa::V = V(undef, N)
    PM::M = M(undef, N,N)
    PS::M = M(undef, N,N)

    Finometrics.project!(X,PM)
    Finometrics.project!(X,Pa)

    Finometrics.projectSlow!(X, PS)
    println("P: ", Pa[1:5])
    println("P (from full matrix): ", PM[1:3,1:3])
    println("P Slow: ", diag(PS)[1:10])
  end
end

function rapidreg(::Type{M}=Matrix{Float64}, ::Type{V}=Vector{Float64};
    iter::Int = 100, N::Int = 10_000, K::Int = 2, testerrors::Bool = true,
    qrtype::Type = M) where {M<:AbstractMatrix, V<:AbstractVector}

  LMtest(M, V, N=N, testerrors=false, K=K, qrtype=qrtype)
end
#LMtest(CuMatrix{Float32}, CuVector{Float32}, N=1_000)
#@time LMtest(CuMatrix{Float32}, CuVector{Float32}, N=500, testerrors=true, K=10)#, qrtype=CuMatrix{Float32})
#CuArrays.allowscalar(false)
#@time LMtest(CuMatrix{Float32}, CuVector{Float32}, N=500, testerrors=true, K=10)#, qrtype=CuMatrix{Float32})

@time rapidreg(CuMatrix{Float32}, CuVector{Float32}, iter=700, N=2_000_000, K=200,
  qrtype=CuMatrix{Float32}) #, qrtype=CuMatrix{Float32})
