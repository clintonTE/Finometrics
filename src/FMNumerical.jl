

#From http://home.online.no/~pjacklam/notes/invnorm/

function FMΦInv(p::Float64)::Float64

    a1::Float64 = -39.6968302866538
    a2::Float64 = 220.946098424521
    a3::Float64 = -275.928510446969
    a4::Float64 = 138.357751867269
    a5::Float64 = -30.6647980661472
    a6::Float64 = 2.50662827745924

    b1::Float64 = -54.4760987982241
    b2::Float64 = 161.585836858041
    b3::Float64 = -155.698979859887
    b4::Float64 = 66.8013118877197
    b5::Float64 = -13.2806815528857

    c1::Float64 = -7.78489400243029E-03
    c2::Float64 = -0.322396458041136
    c3::Float64 = -2.40075827716184
    c4::Float64 = -2.54973253934373
    c5::Float64 = 4.37466414146497
    c6::Float64 = 2.93816398269878

    d1::Float64 = 7.78469570904146E-03
    d2::Float64 = 0.32246712907004
    d3::Float64 = 2.445134137143
    d4::Float64 = 3.75440866190742

    #Rational approximation for lower region.


    if p < 0.02425
        q::Float64 = √(-2.0 * log(p))
        ans::Float64 = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)
    elseif p < 0.97575
        q = (p - 0.5) * (p - 0.5)
        ans = (((((a1 * q + a2) * q + a3) * q + a4) * q + a5) * q + a6) * (p - 0.5) / (((((b1 * q + b2) * q + b3) * q + b4) * q + b5) * q + 1)
    else
        q = √(-2 * log(1 - p))
        ans = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1)
    end

    return ans
end



FMΦ(z::Real)::Real = 0.5+0.5*erf(z*sqrt2Inv)
FMϕ(z::Real)::Real = sqrt2PiInv*exp(-z^2.0/2.0)


#in-place binary search, returns only a boolean value
function inSorted(val::T, B::Vector{T})::Bool where T<:Any


  searched::Bool = false
  low::Int = 1
  high::Int = length(B)
  pivot::T = B[(high+low) ÷ 2]
  found::Bool = false

  #in-place binary sort
  while !found && !searched

    #println("high: $high low: $low mid: $((high+low)÷2) pivot: $pivot")

    if val > pivot
      low = (high+low) ÷ 2 + 1
    elseif val < pivot
      high = (high+low) ÷ 2 - 1
    else
      found = true
    end

    if high-low ≤ 1
      searched = true
      found = val == B[low] || val == B[high]
    end
    pivot = B[(high+low) ÷ 2]
  end

  return found
end

#=function testInSorted(;N::Int = 2000, K::Int = 10_000)

  for i ∈ 1:K
    A::Vector{Int} = rand(1:N, N)
    B::Vector{Int} = rand(1:N, N)

    sort!(B)
    for j ∈ 1:length(A)
      if inSorted(A[j],B) ≠ (A[j] ∈ B)
        error("inSorted:$(inSorted(A[j],B)) while A∈B:$(A[j] ∈ B)\nA: $(A[j])\nB:$B")
      end
    end
  end
end

@time testInSorted()=#

#returns a list of all elements in A
#compare in performance to setdiff. Maintains order unless sort=true
function ANotB(A::Vector{T}, B::Vector{T};
  sort::Bool=false, sorted::Bool=false) where T<:Any

  if !sort && !sorted
    A = deepcopy(A)
    B = deepcopy(B)
  end

  #Force mergesort to deal with sorted cases
  if !sorted
    if !issorted(A)
      sort!(A)
    end
    if !issorted(B)
      sort!(B)
    end
  end

  numA::Int = length(A)
  numB::Int = length(B)

  findings::Vector{Bool} = Vector{Bool}(numA)


  @fastmath @inbounds @simd for i ∈ 1:numA
    findings[i] = length(searchsorted(B,A[i]))>0
  end

  return A[findings]

end

function ANotBTest(N::Int = 1_000_000)
  println("\nStarting test of algorithm N=$N")
  A::Vector{Symbol} = (Symbol).(rand(1:N, N))
  B::Vector{Symbol} = (Symbol).(rand(1:N, N))

  @time setdiff(A,B) #for comparison
  @time ANotB(A,B, sort=true)
end
