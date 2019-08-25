module Finometrics
using Revise

#=if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end=#

#update steps
#cd [current project path]
#] activate .
#] update
#] resolve
#] gc
#] build
#] precompile

##################Dependencies

using  DataFrames, Distributions, StatsBase, GLM, CategoricalArrays,
  Dates, NLopt, ForwardDiff, Formatting, DataStructures, LinearAlgebra

#=function testfunc()
  local f::FormulaTerm

  f = @formula(x ~ y)

  println(f)
end

testfunc()=#

import Base: +, -, ==, >, <, ≥, ≤, isless, isequal
import Base: push!, length, iterate, broadcastable

#######################MACROS###################
#useful macro for conditionally running things in parallel
macro mpar(cond, expr)
  quote
    if $(esc(cond))
        :($(Threads.@threads($expr)))
    else
        :($($expr))
    end
  end
end


######################Methods####################
export FMLM, #Regression methods
  FM2SLS,
  project!,
  getResid!,
  pullModelMatrix,
  getCoeff!,
  getHomoskedΣ!,
  getModWhiteΣ!,
  getWhiteΣ!,
  getNeweyWest!,
  getModWhiteΣSlow,
  getHomoskedΣSlow,
  getWhiteΣSlow,
  getNeweyWestSlow,
  dropNullsFromDF,
  FMExpr,
  FMSym,
  FMQR,
  get1stStage,
  getTerm,
  getR,
  getR²,
  getNeweyWestFunc,
  getClustered!,

  texTable, #IO Mthods
  writeTables2File,
  writeNakedTable,
  array2String,
  num2Str,
  vec2String,
  vbind,

  skewnessStat, #FMstat methods
  kurtosisStat,
  JBStat,
  UnivariateTest,
  autoρ,
  ljungBox,
  test,
  pAutoρ,
  ARMA,
  ARMAFillYHat!,
  ARMALogLike,
  ARMAEstimate,
  generalizeARMA!,
  narrowARMA!,
  winsorize!,
  lagwithin!,
  differencewithin!,

  FMΦInv, #FMNumerical
  FMΦ,

  VAR, #FMVAR
  lagDF!,
  VAREstimate,
  VARPrediction,
  varΣ,
  propagateImpulse,

  FMSpecs, #FMSpecs
  computeFMLMresults!,

  NBool, #Exported types
  NDict,
  NString,
  NSymbol,
  NType,
  NInt,
  NFloat64,
  NDate,

  MString,
  MBool,
  MFloat64,
  MInt,
  MSymbol,
  MDate,
  ∞,

  YearQuarter,
  MYearQuarter,
  YearMonth,
  MYearMonth,
  eom,
  bom,
  eoq,
  boq


##################Custom types

###FMReg Types
const FMExpr = Union{Symbol,Expr,Nothing}
const FMData = Union{Float64,Int,Date,Symbol, CategoricalValue, Missing, Nothing}#,
abstract type FMModel end

###Convenience Types for Export
const NBool = Union{Nothing, Bool}
const NSymbol = Union{Nothing, Symbol}
const NType = Union{Nothing, Type}
const NString = Union{Nothing, String}
const NInt = Union{Int, Nothing}
const NFloat64 = Union{Float64, Nothing}
const NDate = Union{Date, Nothing}

const MBool = Union{Bool, Missing}
const MInt = Union{Int, Missing}
const MFloat64 = Union{Float64, Missing}
const MSymbol = Union{Symbol, Missing}
const MString = Union{String, Missing}
const MDate = Union{Date, Missing}
const ∞ = Inf


###VAR Types
const Structure = Union{Nothing, Vector{Float64}}
const CoefficientErrors = Union{Nothing, Vector{Float64}}
const RegressionSS = Union{Nothing, Float64}
const RegressionN = Union{Nothing, Int}

################Constants

#shared cosntants
const DEFAULT_STAR_LEGEND = "*p<0.1, **p<0.05, ***p<0.01"
const DEBUG_FMMOD = false
const DEFAULT_DECIMALS = 4

#Stat constants
const MAX_TIME = 30.
const F_TOL_ABS = 10. ^ -12.
const F_TOL_REL = 10. ^ -10.
const OPT_BOUND = 2.

#numerical constants
const sqrt2Inv = 1.0/(2.0^.5)
const sqrt2PiInv = 1.0/(2.0*π)^0.5
#  MFloat64, MInt, MBool, MSymbol}

const UPPER_PAD = "%This table was programatically generated using the TexTables package and my own code.\n"
const LOWER_PAD = "\n"

###################Files############

# include("$(pwd())/FMReg.jl")
include("FMReg.jl")
include("FMIO.jl")
include("FMStat.jl")
include("FMVAR.jl")
include("FMNumerical.jl")
include("FMSpecs.jl")
include("FMYearQuarter.jl")
include("FMRegBroken2SLS.jl")

#type is derived from a user-defined type
const MYearQuarter = Union{YearQuarter, Missing}
const MYearMonth = Union{YearMonth, Missing}



end
