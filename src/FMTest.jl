using Revise, Dates, DataFrames, GLM
#NSymbol = Union{Nothing, Symbol}
#NInt = Union{Nothing, Int} #NOTE: consider deleting this later
import Base: +, -, ==, >, <, ≥, ≤, length, isless, isequal
import StatsModels: implicit_intercept


#include("Finometrics.jl")
using Revise
const FMExpr = Union{Symbol,Expr,Nothing}
#  StatsModels.implicit_intercept(::Type{<:Any}) = true
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

function testMM()
  df = DataFrame(x = 1001:2000, y = rand(1000), z=rand(1000),
  f = (i->Symbol(:q,i÷100)).(1:1000))
  #println(df)
  x = getModelMatrix(df, :(x+y))

  println(x)
end

testMM()
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

abstract type YearPeriod end
##################################Year-Month#####################
struct YearMonth <: YearPeriod
  y::UInt16
  m::UInt8

  function YearMonth(y,m)
    (m≥1 && m≤12) || error("Invalid YearMonth $y $m")
    return new(y,m)
  end
end

#creates a year month from yyyymm date style
function YearMonth(i::Integer)::YearMonth
  local y::UInt16
  local m::UInt8

  (i>1000_00 && i<10000_00) || error("Invalid YearMonth $i")
  (y,m) = (i ÷ 100, i % 100)

  return YearMonth(y,m)
end

#fallback in case we need these
Int(ym::YearMonth)::Int = ym.y*100 + ym.m
Base.show(io::IO, ym::YearMonth) = print(io, "$(ym.y)_$(ym.m)")

#quick way to split the year quarter
Tuple{Int,Int}(ym::YearMonth)::Tuple{Int,Int} = (ym.y, ym.m)

#allows for adding quarters
function +(ym::YearMonth, madded::Integer)::YearMonth
  (y₀::Int, m₀::Int) = Tuple{Int,Int}(ym)


  m::Int = (m₀ + madded) % 12
  (m≤0) && (m+=12) #only allow positive months

  y::Int = (madded - (m-m₀)) ÷ 12 + y₀

  @assert m-m₀ + 12*(y-y₀) == madded #this should never fail

  return YearMonth(y,m)
end

#helper methods to build on the above
+(qadded::Int, yq::YearMonth)::YearMonth = ym + qadded
-(ym::YearMonth, madded::Int)::YearMonth = ym + (-1*madded)


#gets the number of months difference
function -(ym1::YearMonth, ym2::YearMonth)::Int

  Δy::Int = ym1.y - ym2.y
  Δm::Int = ym1.m - ym2.m

  return Δy * 12 + Δm
end

#comparison operators (primary sorting functions, then operator symbols)
isless(ym1::YearMonth, ym2::YearMonth) = ym1.y*100+ym1.m < ym2.y*100 + ym2.m#(ym1.y<ym2.y) || ((ym1.y==ym2.y) && (ym1.m<ym2.m))

#<(ym1::YearMonth, ym2::YearMonth) = isless(ym1, ym2)
#>(ym1::YearMonth, ym2::YearMonth) = isless(ym2, ym1)

#≤(ym1::YearMonth, ym2::YearMonth) = isless(ym1, ym2) || isequal(ym1, ym2)
#≥(ym1::YearMonth, ym2::YearMonth) = isless(ym2, ym1) || isequal(ym1, ym2)

#convert back to a date format
bom(ym::YearMonth)::Date = Date(ym.y, ym.m, 1)
eom(ym::YearMonth)::Date = lastdayofmonth(boq(ym))

Base.length(::YearPeriod) = 1
Base.iterate(yp::YearPeriod) = yp
Base.iterate(ym::YearPeriod, s::Integer) = nothing
Base.Broadcast.broadcastable(yp::YearPeriod) = Ref(yp)

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
