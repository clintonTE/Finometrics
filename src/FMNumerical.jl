

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

#vcats like with rbind in r
function vbind(df1::AbstractDataFrame, df2::AbstractDataFrame)::AbstractDataFrame
  #get info on what we are joining
  fields1::Vector{<:DField} = names(df1)
  fields2::Vector{<:DField} = names(df2)
  allFields::Vector{<:DField} = union(fields1, fields2)

  rows1::Int = size(df1,1)
  rows2::Int = size(df2,1)

  for f::DField ∈ setdiff(allFields, fields1)
    origType::Type = eltype(df2[f])

    if !(Missing <: origType) #in this case we need to promote the eltype
      origType = Union{origType, Missing}
      df2[f] = Vector{origType}(df2[f])
    end

    df1[ f] = Vector{origType}(missing, rows1)
  end

  for f::DField ∈ setdiff(allFields, fields2)
    origType = eltype(df1[f])

    if !(Missing <: origType) #in this case we need to promote the eltype
      origType = Union{origType, Missing}
      df1[f] = Vector{origType}(df1[f])
    end

    #println("origType: $origType")
    df2[ f] = Vector{origType}(missing, rows2)
  end

  return [df1; df2]
end

#recursive version that works similar to vcat with the above stipulations
function vbind(dfs::T...)::T where T<:AbstractDataFrame
  if length(dfs) == 1
    return dfs[1]
  else
    pivot::Int = length(dfs) ÷ 2
    return vbind(vbind(dfs[1:pivot]...), vbind(dfs[(pivot+1):end]...))
  end
end
