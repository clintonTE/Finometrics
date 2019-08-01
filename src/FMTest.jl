using DataFrames, Revise
#NSymbol = Union{Nothing, Symbol}
#NInt = Union{Nothing, Int} #NOTE: consider deleting this later
#import Base: +, -, ==, >, <, ≥, ≤


using Finometrics



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
end

#@time for i ∈ 1:200_000
#  testYearQuarter()
#end

function testexisting()
  @assert ismissing(existing(missing))
  @assert existing(1) == 1
  @assert existing(missing,2) == 2
  @assert existing(1,2) == 1
  @assert existing(1,missing) == 1
  @assert existing(1,2,3) == 1
  @assert existing(missing,2,3) == 2
  @assert existing(1,missing,3) == 1
  @assert existing(1,2,missing) == 1
  @assert existing(missing,2,missing) == 2
  @assert existing(missing,missing,3) == 3
  @assert existing([missings(100);101]...) == 101
  @assert existing([1; missings(100)]...) == 1
  @assert existing([missings(50); 51; missings(50)]...) == 51
  @assert existing(collect(1:101)...) == 1
end

function testexisting2()
  @assert ismissing(coalesce(missing))
  @assert coalesce(1) == 1
  @assert coalesce(missing,2) == 2
  @assert coalesce(1,2) == 1
  @assert coalesce(1,missing) == 1
  @assert coalesce(1,2,3) == 1
  @assert coalesce(missing,2,3) == 2
  @assert coalesce(1,missing,3) == 1
  @assert coalesce(1,2,missing) == 1
  @assert coalesce(missing,2,missing) == 2
  @assert coalesce(missing,missing,3) == 3
  @assert coalesce([missings(100);101]...) == 101
  @assert coalesce([1; missings(100)]...) == 1
  @assert coalesce([missings(50); 51; missings(50)]...) == 51
  @assert coalesce(collect(1:101)...) == 1
end

@time testexisting()

#NOTE: Originally tried to create a primitive Quarter type,
#keep as reference
#=
primitive type YearQuarter <: AbstractFloat 64 end


function YearQuarter(y::Int, q::Int)
  #default effectively means that a provided year is treated as the fourth quarter of that year
  q == 1 || q==2 || q==3 || q==4 || error("Invalid YearQuarter $f")

  return reinterpret(YearQuarter, y + q/10)
end


function YearQuarter(f::Real)::YearQuarter
  y::Int = round(f)
  q::Int = (f*10-y*10)
  (q==0) && (q=4) #if a year is passed, assume it is the end of the year (q4)
  q == 1 || q==2 || q==3 || q==4 || error("Invalid YearQuarter $f")
  return reinterpret(YearQuarter, f)
end

#fallback in case we need these
Float64(yq::YearQuarter) = reinterpret(Float64, yq)
Base.show(io::IO, yq::YearQuarter) = print(io, Float64(yq))

#quick way to split the year quarter
function Tuple{Int,Int}(yq::YearQuarter)::Tuple{Int,Int}

  f::Float64 = Float64(yq)
  y::Int = Int(round(f))
  q::Int = f*10 - y*10

  return (y,q)
end

#allows for adding quarters
function +(yq::YearQuarter, qadded::Int)::YearQuarter
  (y::Int, q::Int) = yq

  y += (qadded) ÷ 4
  rawq::Int = qadded - yadded * 4 + q #should be ∈ -2:7

  #this handles the edge cases for the quarters
  if (rawq > -3) && (rawq ≤ 0)
    y -= 1
    q = rawq + 4
  elseif (1≤rawq) && (rawq≤4)
    q = rawq
  elseif (5≤rawq) && (rawq≤7)
    y += 1
    q = rawq - 4
  else
    @assert false
  end

  return reinterpret(YearQuarter, y + q/10)
end

#helper methods to build on the above
+(qadded::Int, yq::YearQuarter)::YearQuarter = yq + qadded
-(yq::YearQuarter, qadded::Int)::YearQuarter = yq + (-1*qadded)

#gets the number of quarters difference
function -(yq1::YearQuarter, yq2::YearQuarter)::Int
  (y1::Int, q1::Int) = yq1
  (y2::Int, q2::Int) = yq2

  Δy::Int = y1 - y2
  Δq::Int = q1 - q2

  return Δy * 4 + Δq
end

function testYearQuarter()

  #test the first constructor
  yq1::YearQuarter = YearQuarter(1985,3)
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

  println(yq1)
end

testYearQuarter()=#
