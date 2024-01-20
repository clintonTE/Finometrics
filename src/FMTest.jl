#NOTE: Uncomment below bloc for stand alone testing
#WARNING WARNING WARNING don't for get to comment include("test.jl") in Finometics.jl
using Revise

#=include("Finometrics.jl")
using Distributions, LinearAlgebra, #=CUDA, =#DataFrames,
 Dates, DataFrames, GLM, Random
Random.seed!(11)
import Base: +, -, ==, >, <, ≥, ≤, length, isless, isequal
import StatsModels: implicit_intercept

const MFloat64 = Union{Float64, Missing}
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
end=#


#NOTE: End stand-alone comment block

function testMM(df::DataFrame, loops=1)
  #println(df)
  m = Finometrics.generateX(Matrix{Float64}, df, :(x+y+f))

  (loops==1) && display(m)

  return m
end

function testMM(loops::Int; N::Int=1000, G::Int=50)
  local tot::Float64 = 0.0
  @info "testing MM. Failure will result in an error."

  df = DataFrame(x = 1001:(1000+N), y = rand(N), z=rand(N),
    f = (i->Symbol(:q,i÷(N ÷ G))).(1:N))

  @time for i ∈ 1:loops
    m=testMM(df, loops)
    (sum(df.x .≠ m[:,2]) == 0) || error("Invalid x column! df.x:$(df.x)\n###\nm[:,2]:$(m[:,2])")
    (sum(df.y .≠ m[:,3]) == 0) || error("Invalid y column!")
    (sum(1 .≠ m[:,1]) == 0) || error("Invalid intercept column!")
    (size(m,1) == N) || error("Invalid row dimension!")
    (size(m,2) == length(unique(df.f))+2) || error(
      "Invalid column dimension! Effects=$(length(unique(df.f))) size=$(size(m))")
  end

  println(tot)
end


#use this to test the lag functions
function testlaganddifference(N = 10_000)
  rowindex = collect(1:N)
  periodcol = deepcopy(rowindex)

  @info "testing lag and difference. Failure will result in an error."


  #make a test frame
  df = DataFrame(rowindex=rowindex, period = periodcol, group1 = periodcol .÷ 17,
    group2 = periodcol .÷ 29, target1 = rand(N), target2 = rand(N))

  #sort!(df, [:group1, :period, :group2])
  #CASE 1: create lagged field Lgroup1
  Finometrics.lagwithin!(df, :target1, :group1, :period)

  #CASE 2: create lagged fields L2target1 and L2target2
  sort!(df, [:group1, :group2])
  Finometrics.lagwithin!(df, [:target1, :target2], [:group1, :group2], :period, lags=2)

  #CASE 3: difference :Ltarget1 and create DLtarget1 and LLtarget1
  Finometrics.differencewithin!(df, :Ltarget1, :group1, :period, deletelag=false,
    laggedname=:LLtarget1e, differencedname=:DLtarget1e)

  #CASE 4: difference [:L2target1, :L2target2] creating
  # [:DL2target1e, :DL2target2e], [:LL2target1e, :LL2target2e]
  Finometrics.differencewithin!(df, [:L2target1, :L2target2], [:group1, :group2], :period,
    deletelag=false, differencednames=[:DL2target1e, :DL2target2e],
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


function testYearQuarter()

  @info "testing year quarter. Failure will result in an error."
  #test the first constructor
  yq1::Finometrics.YearQuarter = Finometrics.YearQuarter(1998,2)
  try
    YearQuarter(1985,7)
    @assert false
  catch err
    (typeof(err)<:AssertionError) && error("Error catching invalid values in first constructor")
  end

  #test the second constructor
  yq2::Finometrics.YearQuarter = Finometrics.YearQuarter(2000.4)
  try
    Finometrics.YearQuarter(2000.05)
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

function testYearMonth()
  @info "testing year month. Failure will result in an error."
  #test the first constructor
  ym1::Finometrics.YearMonth = Finometrics.YearMonth(1998,2)
  try
    Finometrics.YearMonth(1985,14)
    @assert false
  catch err
    (typeof(err)<:AssertionError) && error("Error catching invalid values in first constructor")
  end

  #test the second constructor
  ym2::Finometrics.YearMonth = Finometrics.YearMonth(2008_04)
  try
    Finometrics.YearMonth(2000_14)
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



function LMtest(::Type{M}=Matrix{Float64}, ::Type{V}=Vector{Float64};
    N::Int = 200, K::Int = 2, testerrors::Bool = true,
    qrtype::Type = M, iter::Int = 1, testprimarywithin::Bool = false, runslow::Bool = true) where {
    M<:AbstractMatrix, V<:AbstractVector}

  @warn "Testing LM. Most tests are automated at this point but
    some tests, like cluster equivelency, are manual. Generally though, for a quick look
    a lack of error should indicate nothing horrific is going on."
  local lin::Finometrics.FMLM

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
  for i ∈ 1:K
    X[:,i] .= V(rand(Uniform(),N))
  end

  ε = ((s2::Float64)->rand(Normal(0.0,s2^0.5))).(σ2)
  Y =  X*β .+ ε

  df::DataFrame = DataFrame(idx = 1:N)
  xnames = (i->Symbol(:X,i)).(1:K)

  for (i,x) ∈ enumerate(eachcol(X))
    df[!, xnames[i]] = x
  end
  df.Y = Y
  df.G = (i->Symbol(:G,i)).(rand(1:5,N))
  df.C1 = (i->Symbol(:C,i)).((j->j÷10).(1:N))
  df.C2 = (i->Symbol(:C,i)).((j->j%10).(1:N))

  #get the linear model
  xspec = Meta.parse(join((string).(xnames),"+"))

  println("minimal test")
  Finometrics.FMLM(df, xspec, :Y, containsmissings=false)

  print("Time running primary regression...")
  @time for i ∈ 1:iter
    lin = Finometrics.FMLM(df, xspec, :Y, M, V, #withinsym=:C1,
      clustersyms=:C1, qrtype=qrtype, checkwithin=testprimarywithin, containsmissings=false)
  end

  if testerrors
    #get the homoskedastic SEs
    ΣHomosked::Matrix{Float64} = Finometrics.homoskedasticΣ!(lin)
    runslow && (ΣHomoskedSlow::Matrix{Float64} = Finometrics.homoskedasticΣslow(lin))


    #print the coefficients
    println("\nCoefficients: ",lin.β)
    modelcoef = ((lin.X' * lin.X)\I) * lin.X' * lin.Y
    println("Model coeff: ", modelcoef)
    @assert modelcoef ≈ lin.β
    A::Matrix{Float64} = [ones(N) df.X1 df.X2 df.X3 df.X4 df.X5]
    manualcoef = ((A' * A)\I) * A' * Y
    @assert manualcoef ≈ lin.β
    println("Manual coef:", manualcoef)
    @info "passed coefficient test"

    xnamesfixed = [xnames; :G]
    xspecfixed = Meta.parse(join((string).(xnamesfixed),"+"))
    linfixed = Finometrics.FMLM(df, xspecfixed, :Y, M, V, #withinsym=:C1,
      clustersyms=[:C1, :C2], qrtype=qrtype, checkwithin=testprimarywithin, containsmissings=false)
    xspecwithin = Meta.parse("$(join((string).(xnames),"+")) + 0")
    linwithin::Finometrics.FMLM = Finometrics.FMLM(df, xspecwithin, :Y, M, V, withinsym=:G,
      clustersyms=[:C1, :C2], checkwithin=true, qrtype=qrtype)
    ΣHomoskedfixed::Matrix{Float64} = Finometrics.homoskedasticΣ!(linfixed)
    ΣHomoskedwithin::Matrix{Float64} = Finometrics.homoskedasticΣ!(linwithin)

    println("Coefficients (fixed G): ",linfixed.β)
    println("Coefficients (w/in G): ",linwithin.β)
    @assert linwithin.β ≈ linfixed.β[2:(K+1)]
    @info "passed coefficients fixed effect vs within test"

    homoskederr = diag(ΣHomosked).^.5
    homoskederrslow = diag(ΣHomoskedSlow).^.5
    println("\nHomoskedastic Errors: ", homoskederr)
    runslow && println("Check: ", homoskederrslow)
    runslow && (@assert homoskederr ≈ homoskederrslow)
    runslow && @info "passed homosked test"

    homoskederrfixed = diag(ΣHomoskedfixed)[2:(1+K)].^.5
    homoskederrwithin = diag(ΣHomoskedwithin).^.5
    println("Homoskedastic Errors (fixed): ", homoskederrfixed)
    println("Homoskedastic Errors (w/in): ", homoskederrwithin)
    @assert homoskederrfixed ≈ homoskederrwithin
    @info "passed homosked fixed effects vs within test"

    @info "check for warnings if intercept + within"
    linwithin = Finometrics.FMLM(df, xspec, :Y, M, V, withinsym=:G,
      clustersyms=[:C1, :C2], checkwithin=true, qrtype=qrtype)

    #get the modified white SEs
    ΣMWhite::Matrix{Float64} = similar(ΣHomosked)
    Finometrics.modifiedwhiteΣ!(lin, ΣMWhite)
    runslow && (ΣMWhiteSlow::Matrix{Float64} = Finometrics.modifiedwhiteΣslow(lin))

    #print the coefficients
    println("\nModified White Errors: ",diag(ΣMWhite).^.5)
    runslow && println("Check: ", diag(ΣMWhiteSlow).^.5)
    runslow && (@assert diag(ΣMWhite).^.5 ≈ diag(ΣMWhiteSlow).^.5)
    runslow && @info "passed MWhite test"

    #get the white SEs
    ΣWhite::Matrix{Float64} = similar(ΣHomosked)
    Finometrics.whiteΣ!(lin, ΣWhite)
    runslow && (ΣWhiteSlow::Matrix{Float64} = Finometrics.whiteΣslow(lin))

    #print the coefficients
    println("\nWhite Errors: ",diag(ΣWhite).^.5)
    runslow && println("Check: ", diag(ΣWhiteSlow).^.5)
    runslow && (@assert diag(ΣWhite).^.5 ≈ diag(ΣWhiteSlow).^.5)
    runslow && @info "passed white error test"

    

    #get the modified white SEs
    Σclustered::Matrix{Float64} = similar(ΣHomosked)
    Finometrics.clusteredΣ!(lin, Σclustered)
    ΣclusteredChk =Finometrics.clusteredΣ!(linwithin, Σclustered) #might need a better check here
    println("\nClustered Errors: ",diag(Σclustered).^.5)
    println("Check: ", diag(ΣclusteredChk ).^.5)
    @assert diag(Σclustered).^.5 ≈ diag(ΣclusteredChk ).^.5
    @info "passed cluster error test"

    linalt2 = Finometrics.FMLM(df, xspecwithin, :Y, M, V, withinsym=:G,
      clustersyms=[:C1,:C2], qrtype=qrtype, checkwithin=testprimarywithin, containsmissings=false)
    Σwithinclustered::M = M(undef, linalt2.K, linalt2.K)
    Finometrics.clusteredΣ!(linalt2, Σwithinclustered, clusters = [linalt2.clusters[1]], testequivelance=true)
    Finometrics.clusteredΣ!(linalt2, Σwithinclustered, testequivelance=true)
    @info "passsed cluster method equivelance test"
    

    #get the nw SEs
    ΣNW::Matrix{Float64} = similar(ΣHomosked)
    Finometrics.neweywestΣ!(lin, 3, ΣNW)
    runslow && (ΣNWSlow::Matrix{Float64} = Finometrics.neweywestΣslow(lin, 3))
 
    #print the coefficients
    println("\nNW Errors: ",diag(ΣNW).^.5)
    runslow && println("Check: ", diag(ΣNWSlow).^.5)
    runslow && @assert diag(ΣNW).^.5 ≈ diag(ΣNWSlow).^.5
    runslow && @info "passed newey west test"
      

    #get the nw SEs
    ΣNWpanel::Matrix{Float64} = similar(ΣHomosked)
    Finometrics.neweywestpanelΣ!(lin, 3, ΣNWpanel)
    #throw("got here2.5") 
    runslow && (ΣNWpanelslow::Matrix{Float64} = Finometrics.neweywestpanelΣslow!(lin, 3))
    #throw("got here3") 
    
    #print the coefficients
    println("\nNW panel Errors: ",diag(ΣNWpanel).^.5)
    runslow && println("Check: ", diag(ΣNWpanelslow).^.5)
    runslow && diag(ΣNWpanel).^.5 ≈ diag(ΣNWpanelslow).^.5
    runslow && @info "passed newey west panel test"
    #throw("got here3") 

    SST::Float64 = sum((lin.Y .- mean(lin.Y)).^2)
    SSE::Float64 = sum(lin.ε.^2)

    println("\nR²: $(Finometrics.R²(lin)) R²-derived: $(1 - SSE/SST)")
    @assert Finometrics.R²(lin) ≈ 1 - SSE/SST
    @info "passed R² test"
    
    #test the project routines
    Pa::V = V(undef, N)
    Finometrics.project!(X,Pa)
    println("\nP: ", Pa[1:5])
    if N≤10_000 #only run on small samples or we will run out of memory
      PM::M = M(undef, N,N)
      runslow && (PS::M = M(undef, N,N))

      Finometrics.project!(X,PM)
      runslow && Finometrics.projectslow!(X, PS)
      println("P (from full matrix): ", PM[1:3,1:3])
      @assert diag(PM) ≈ Pa
      runslow && println("P Slow: ", diag(PS)[1:10])
      runslow && @assert diag(PS) ≈ Pa
      @info "projection tests passed"
    end
  end
end

function rapidreg(::Type{M}=Matrix{Float64}, ::Type{V}=Vector{Float64};
    iter::Int = 100, N::Int = 1_000, K::Int = 2, testerrors::Bool = true,
    qrtype::Type = M, testprimarywithin=false) where {M<:AbstractMatrix, V<:AbstractVector}

  LMtest(M, V, N=N, testerrors=false, testprimarywithin=testprimarywithin, K=K, qrtype=qrtype, iter=iter)
end

function testlagwithin2(N=100_000)
  @info "testing lagwithin2. Failure will result in an error."
  #make the test data
  df::DataFrame = DataFrame(
    val1 = Vector{Union{Float64, Missing}}(undef, N),
    val2 = Vector{Union{Float64, Missing}}(undef, N),
    group = rand(collect(1:1000),N),
    date = rand(collect(Date(1985,11,11):Day(1):Date(2011,11,11)), N))

  for i ∈ 1:N, s ∈ [:val1, :val2]
    if rand() < 0.99
      df[i, s] = rand()
    else
      df[i, s] = missing
    end
  end

  maxnotstale = Day(100)

  #execute! (will also sort)
  try
    Finometrics.lagwithin2!(df, [:val1, :val2],  :group, date=:date) #minimalist version, unsorted
    @assert false
  catch err
    typeof(err) <: AssertionError && error("lagwithin2 should have failed due to lack of sorting")
    sort!(df, [:group, :date])
  end
  dfnostale = deepcopy(df)

  Finometrics.lagwithin2!(dfnostale, [:date, :val1, :val2],  :group, date=:date) #minimalist version, sorted
  Finometrics.leadwithin2!(dfnostale, [:date, :val1, :val2],  :group, date=:date) #minimalist version, sorted
  Finometrics.lagwithin2!(df, [:date, :val1, :val2],
    :group, date=:date, maxnotstale = maxnotstale)
  Finometrics.leadwithin2!(df, [:date, :val1, :val2],
    :group, date=:date, maxnotstale = maxnotstale)

  @assert sum(ismissing.(dfnostale.Ldate)) < sum(ismissing.(df.Ldate))
  @assert sum(ismissing.(dfnostale.Lval1)) < sum(ismissing.(df.Lval1))
  @assert sum(ismissing.(dfnostale.Lval2)) < sum(ismissing.(df.Lval2))
  @assert sum(ismissing.(dfnostale.Ndate)) < sum(ismissing.(df.Ndate))
  @assert sum(ismissing.(dfnostale.Nval1)) < sum(ismissing.(df.Nval1))
  @assert sum(ismissing.(dfnostale.Nval2)) < sum(ismissing.(df.Nval2))

  #Finometrics.leadwithin2!(df, [:val1, :val2],  :group, date=:date) #minimalist version, sorted
  #Finometrics.lagwithin2!(df, [:val1, :val2],  :group, date=:date) #minimalist version, sorted

  #test the no dataframe version
  df.vecLval1 = Finometrics.lagwithin2sorted(df.val1,df.group)
  @assert all(dfnostale.Lval1 .=== df.vecLval1)

  Finometrics.differencewithin2!(df, [:val1, :val2],
    :group, date=:date, maxnotstale = maxnotstale)

  Finometrics.leadwithin2!(df, [:date, :val1, :val2],
    :group, date=:date, maxnotstale = maxnotstale)

  Finometrics.lagwithin2!(df, [:Ndate, :Nval1, :Nval2],
    :group, date=:date, maxnotstale = maxnotstale)

  Finometrics.leadwithin2!(df, [:Ldate, :Lval1, :Lval2],
    :group, date=:date, maxnotstale = maxnotstale)
  @assert (sum(df.LNval1 .=== df.val1) > 0) && (sum(df.LNval2 .=== df.val2) > 0)

  df.TLdate = similar(df.Ldate)
  df.TLval1 = similar(df.Lval1)
  df.TLval2 = similar(df.Lval2)
  #make the test lags
  for sdf ∈ groupby(df, :group)
    for j ∈ 2:size(sdf,1)
      if sdf[j-1,:date] ≥ sdf[j,:date] - maxnotstale
        sdf[j,:TLdate] = sdf[j-1,:date]
        sdf[j,:TLval1] = sdf[j-1,:val1]
        sdf[j,:TLval2] = sdf[j-1,:val2]
      end
    end
  end

  #compute the period
  df.idx = 1:N |> collect

  df.daysfromlag = ((d,Ld)->((d===missing) | (Ld === missing) ? missing : d-Ld)
    ).(df.date, [missing; df.date[1:(end-1)]])
  df.Lgroup = [missing; df.group[1:(end-1)]]
  df.daysfromlag[df.Lgroup .!== df.group] .= missing

  df.daystolead = ((d,Nd)->((d===missing) | (Nd === missing) ? missing : Nd-d)
    ).(df.date, [df.date[2:end]; missing])
  df.Ngroup = [df.group[2:end]; missing]
  df.daystolead[df.Ngroup .!== df.group] .= missing

  @assert all(df.Ldate .=== df.TLdate )

  @assert all(df.Lval1 .=== df.TLval1 )
  @assert all((df.val1 .- df.TLval1) .=== df.Dval1)

  @assert all(df.Lval2 .=== df.TLval2 )
  @assert all((df.val2 .- df.TLval2) .=== df.Dval2)

  df.TNdate = similar(df.Ndate)
  df.TNval1 = similar(df.Nval1)
  df.TNval2 = similar(df.Nval2)
  #make the test lags
  for sdf ∈ groupby(df, :group)
    for j ∈ 1:(size(sdf,1)-1)
      if sdf[j,:date] ≥ sdf[j+1,:date] - maxnotstale
        sdf[j,:TNdate] = sdf[j+1,:date]
        sdf[j,:TNval1] = sdf[j+1,:val1]
        sdf[j,:TNval2] = sdf[j+1,:val2]
      end
    end
  end

  @assert all(df.Ndate .=== df.TNdate )
  @assert all(df.Nval1 .=== df.TNval1 )
  @assert all(df.Nval2 .=== df.TNval2 )

  ####NOTE- commentary for testing if leading a lag is equivent to the original
  #this basically says that leading and then lagging the value produces the original value
  #that is, either val = LNval or a missing value is resulting from one of two places:
  #1) the original val is missing 2) its an endpoint, so the Lval1 is missing
  #(we don't need to worry about the Nval missing since the NL value at t draws from the
  #N value at t-1 which draws from the t value- but this notably implies NLVal1=Val1=missing so
  #we don't need an extra case)
  # the NLval=val is handled analagously, where the original value is missing or its
  #an end point so the Nval1 is missing. We don't need to worry about the Lval missing since the NL
  #value draws from the N value of L which draws from the val, which would imply LNVal1=Val1=missing

  #println(df[1:200,[:idx, :date, :group, :val1, :Lval1, :Nval1, :LNval1, :daysfromlag, :daystolead]])
  @assert all(((df.LNdate .=== missing) .& (df.Ldate .=== missing))
   .| (df.LNdate .=== df.date))
  @assert all(((df.LNval1 .=== missing) .& (df.Lval1 .=== missing))
   .| (df.LNval1 .=== df.val1))
   @assert all(((df.LNval2 .=== missing) .& (df.Lval2 .=== missing))
    .| (df.LNval2 .=== df.val2))

  @assert all(((df.NLdate .=== missing) .& (df.Ndate .=== missing))
   .| (df.NLdate .=== df.date))
  @assert all(((df.NLval1 .=== missing) .& (df.Nval1 .=== missing))
   .| (df.NLval1 .=== df.val1))
   @assert all(((df.NLval2 .=== missing) .& (df.Nval2 .=== missing))
    .| (df.NLval2 .=== df.val2))

end


function testwinsorizequantile(i=100_000; tol = 10^-8)
  v::Vector{MFloat64} = rand(i)
  v = (f->rand()<0.05 ? missing : f).(v)

  @inline completev(v) = view(v, (!ismissing).(v))

  println("test: winsorizequantile(v, .8, twotailed=false)")
  wv = completev(Finometrics.winsorizequantile(v, .8, twosided=false))
  println("results (exp=~0.0,0.8): min: $(minimum(wv)), max: $(maximum(wv))")

  println("test: winsorizequantile(v, .3, twotailed=false)")
  wv = completev(Finometrics.winsorizequantile(v, .3, twosided=false))
  println("results (exp=~0.3,1.0): min: $(minimum(wv)), max: $(maximum(wv))")

  println("test: winsorizequantile(v, .2, twotailed=true)")
  wv = completev(Finometrics.winsorizequantile(v, .2, twosided=true))
  println("results (exp=~0.2,0.8) min: $(minimum(wv)), max: $(maximum(wv))")

  println("autotest: winsorizequantile(v, .8, twotailed=true)")
  wv2 = completev(Finometrics.winsorizequantile(v, .2, twosided=true))
  (sum((abs).(wv2 .- wv) .> tol) > 0) && error("symmetric two-tailed winsorization test failed")

  return nothing
end

function runbasictests()
  @info "Some basic tests. Incomplete, but better than nothing until I get around to making something better"

  testMM(100, N=1000, G=10)
  @time LMtest(Matrix{Float64}, Vector{Float64},
    N=2_000, testerrors=true, K=5, testprimarywithin=true, runslow=true)#, qrtype=CuMatrix{Float32})

  testYearQuarter()
  testYearMonth()
  testlaganddifference()
  testlagwithin2()
  testwinsorizequantile()
end


#runbasictests()
