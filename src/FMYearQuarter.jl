struct YearQuarter
  y::UInt16
  q::UInt8

  function YearQuarter(y,q)
    q == 1 || q==2 || q==3 || q==4 || error("Invalid YearQuarter $f")
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
function pold(yq::YearQuarter, qadded::Int)::YearQuarter
  (y::Int, q::Int) = Tuple{Int,Int}(yq)

  yadded::Int = (qadded) ÷ 4
  y::Int += yadded
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

  return YearQuarter(y,q)
end

#allows for adding quarters
function +(yq::YearQuarter, qadded::Int)::YearQuarter
  (y₀::Int, q₀::Int) = Tuple{Int,Int}(yq)


  q::Int = (q₀ + qadded) % 4
  (q≤0) && (q+=4) #only allow 1,2,3,4

  y::Int = (qadded - (q-q₀)) ÷ 4 + y₀

#  @assert q-q₀ + 4*(y-y₀) == qadded #this should never fail

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
isequal(yq1::YearQuarter, yq2::YearQuarter) = (yq1.y==yq2.y) && (yq1.q==yq2.q)
isless(yq1::YearQuarter, yq2::YearQuarter) = (yq1.y<yq2.y) || ((yq1.y==yq2.y) && (yq1.q<yq2.q))

==(yq1::YearQuarter, yq2::YearQuarter) = isequal(yq1, yq2)
<(yq1::YearQuarter, yq2::YearQuarter) = isless(yq1, yq2)
>(yq1::YearQuarter, yq2::YearQuarter) = isless(yq2, yq1)

≤(yq1::YearQuarter, yq2::YearQuarter) = isless(yq1, yq2) || isequal(yq1, yq2)
≥(yq1::YearQuarter, yq2::YearQuarter) = isless(yq2, yq1) || isequal(yq1, yq2)
