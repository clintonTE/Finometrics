

abstract type YearPeriod end

struct YearQuarter <: YearPeriod
  y::UInt16
  q::UInt8

  function YearQuarter(y,q)
    q == 1 || q==2 || q==3 || q==4 || error("Invalid YearQuarter $y $q")
    return new(y,q)
  end
end

#creates a year quarter from the yyyy.q format or from yyyy (assuming eoy)
#unpredicatable results if the decimal is not 0 (year only case), 1,2,3, or 4
function YearQuarter(f::Real)::YearQuarter
  ix10::Int = Int(f*10)
  y::Int = ix10 ÷ 10
  q::Int = ix10 - y*10
  (q==0) && (q=4) #if a year is passed, assume it is the end of the year (q4)

  return YearQuarter(y,q)
end

function YearQuarter(d::Date)
  y::UInt16 = year(d)
  q::UInt8 = quarterofyear(d)

  return YearQuarter(y,q)
end

#fallback in case we need these
Float64(yq::YearQuarter) = yq.y + yq.q/10
Base.show(io::IO, yq::YearQuarter) = print(io, Float64(yq))

#quick way to split the year quarter
Tuple{Int,Int}(yq::YearQuarter)::Tuple{Int,Int} = (yq.y, yq.q)

#allows for adding quarters
function +(yq::YearQuarter, qadded::Int)::YearQuarter
  (y₀::Int, q₀::Int) = Tuple{Int,Int}(yq)


  q::Int = (q₀ + qadded) % 4
  (q≤0) && (q+=4) #only allow 1,2,3,4

  y::Int = (qadded - (q-q₀)) ÷ 4 + y₀

  @assert q-q₀ + 4*(y-y₀) == qadded #this should never fail

  return YearQuarter(y,q)
end

#helper methods to build on the above
+(qadded::Int, yq::YearQuarter)::YearQuarter = yq + qadded
-(yq::YearQuarter, qadded::Int)::YearQuarter = yq + (-1*qadded)

#gets the number of quarters difference
function -(yq1::YearQuarter, yq2::YearQuarter)::Int

  Δy::Int = yq1.y - yq2.y
  Δq::Int = yq1.q - yq2.q

  return Δy * 4 + Δq
end

#comparison operators (primary sorting functions, then operator symbols)
isless(yq1::YearQuarter, yq2::YearQuarter) = yq1.y*10+yq1.q < yq2.y*10 + yq2.q

boq(yq::YearQuarter)::Date = Date(yq.y, (yq.q-1)*3,1)
eoq(yq::YearQuarter)::Date = lastdayofquarter(boq(yq))


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
isless(ym1::YearMonth, ym2::YearMonth) = ym1.y*100+ym1.m < ym2.y*100 + ym2.m

#convert back to a date format
bom(ym::YearMonth)::Date = Date(ym.y, ym.m, 1)
eom(ym::YearMonth)::Date = lastdayofmonth(bom(ym))

Base.length(::YearPeriod) = 1
Base.iterate(yp::YearPeriod) = yp
Base.iterate(ym::YearPeriod, s::Integer) = nothing
Base.Broadcast.broadcastable(yp::YearPeriod) = Ref(yp)
