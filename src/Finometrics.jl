module Finometrics
using Revise


##################Dependencies

using  DataFrames, Distributions, StatsBase, GLM, CategoricalArrays,  Dates, LinearAlgebra, StatsModels, NLopt

try
  using CUDA
  CUDA.allowscalar(false) #default
catch
  println("Note: CuArrays not installed")
  CUDA = Nothing #allows checks for CuArrays to work (should all return false)
end

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
#BROKEN
#=macro mpar(cond, expr)
  quote
    if $(esc(cond))
        :($(Threads.@threads($expr)))
    else
        :($($expr))
    end
  end
end=#


######################Methods####################
export FMLM, #Regression methods
  project!, #regression methods (won't work)
  resid!,
  coeff!,
  homoskedasticΣ!,
  modifiedwhiteΣ!,
  whiteΣ!,
  neweywestΣ!,
  homoskedasticΣslow,
  modifiedwhiteΣslow,
  whiteΣslow,
  neweywestΣslow,
  FMExpr,
  FMSym,
  FMQR,
  βterm, #odd naming is due to conflicts with other packages
  R,
  R²,
  neweywestΣfunc,
  neweywestΣfuncslow,
  neweywestpanelΣ!,
  neweywestpanelΣslow!,
  neweywestpanelΣfunc,
  neweywestpanelΣslowfunc,
  clusteredΣ!,



  textable, #IO Mthods
  writetextable,
  array2string,
  num2str,

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
  winsorizequantile,
  winsorizelevel,

  lagwithin!, #FMLagWithin
  differencewithin!,
  differencewithin2!,
  lagwithin2!,
  lagwithin2sorted!,
  lagwithin2sorted,
  leadwithin2!,

  FMΦInv, #FMNumerical
  FMΦ,
  inSorted,
  vbind,

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
  MAString,
  MBool,
  MFloat64,
  MInt,
  MSymbol,
  MDate,
  ∞,

  CU32,
  CU64,

  DField,
  MDField,
  NDField,

  #yearmonth types
  YearQuarter,
  MYearQuarter,
  YearMonth,
  MYearMonth,
  eom,
  bom,
  eoq,
  boq,

  ##WARNING DEPRECIATED methods #WARNING
  texTable,
  writeTexTable,
  array2String,
  num2Str,
  vec2String,
  getCoeff!,
  getHomoskedΣ!,
  getModWhiteΣ!,
  getWhiteΣ!,
  getNeweyWest!,
  getModWhiteΣSlow,
  getHomoskedΣSlow,
  getWhiteΣSlow,
  getNeweyWestSlow,
  getTerm,
  getR,
  getR²,
  getNeweyWestFunc,
  getClustered!,
  winsorize!

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
const MAString = Union{AbstractString, Missing}
const MDate = Union{Date, Missing}
const ∞ = Inf

const DField = Union{String, Symbol}
const MDField = Union{String, Symbol, Missing}
const NDField = Union{String, Symbol, Nothing}


###VAR Types
const Structure = Union{Nothing, Vector{Float64}}
const CoefficientErrors = Union{Nothing, Vector{Float64}}
const RegressionSS = Union{Nothing, Float64}
const RegressionN = Union{Nothing, Int}

################Constants
const spinlock = Threads.SpinLock()

#shared cosntants
const DEFAULT_STAR_LEGEND = raw"$\text{*}p<0.1$, $\text{**}p<0.05$, $\text{***}p<0.01$"
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
include("FMCov.jl")
include("FMIO.jl")
include("FMStat.jl")
include("FMLagWithin.jl")
include("FMVAR.jl")
include("FMNumerical.jl")
include("FMSpecs.jl")
include("FMYearQuarter.jl")
include("FMTest.jl") #comment this out when testing
#include("FMRegBroken2SLS.jl")

#type is derived from a user-defined type
const MYearQuarter = Union{YearQuarter, Missing}
const MYearMonth = Union{YearMonth, Missing}



end
