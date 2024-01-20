######################IO################

vec2String(args...; keyargs...) = error("Depreciated- replace with join")
texTable(args...; keyargs...) = error("Depreciated- replace with textable")
writeTexTable(args...; keyargs...) = error("Depreciated- replace with writetextable")
array2String(args...; keyargs...) = error("Depreciated- replace with array2string")

#functions for converting Floats to Strings
function num2str(x::T, decimals::Int=DEFAULT_DECIMALS;
    scalefactor::Real = 1.0,
    scalefactoronhurdle::Real = 1.0,
    Ints::Bool = false, scalehurdle::U=nothing, texequation::Bool = false)::String where {
      T<:Union{Real,Missing,Nothing}, U<:Union{Nothing,Real}}

  local outstr::String

  x*=scalefactor

  #this controls for both missing and nothing values
  if ismissing(something(x,missing))
    outstr = ""
  elseif scalehurdle ≠ nothing && x≥scalehurdle
    outstr = "$(Int(round(x*scalefactoronhurdle)))"
  elseif decimals == 0 || (Ints && round(x) == x)
    outstr = "$(Int(round(x)))"
  else
    outstr = format("{:.$(decimals)f}",x)
  end

  outstr = texequation ? "\$$outstr\$" : outstr

  return outstr
end

#depreciated functions
#array2String
  #=takes content of an array, formats each cell into a string, appends a
  prefix and suffix, and sends back an array of strings with the same
  dimension as the input
  IN: A general array, optionally a prefix, suffix, scaling factor, and
   set number of decimal digits
  OUT: An array of strings with the same dimensions as the input=#
function array2string(m::Array{Float64}; prefix::String = "",
  suffix::String = "", scaling::Float64 = 1.0,
  decimaldigits::Int = 2)::Array{String}

  #broadcast the prefix, value in m, and suffix and concatenate respectively
  stringmat::Array{String} =
      ((x::Float64)->(prefix*"$(round(x,decimaldigits))"*suffix)).(m)

  return stringmat
end

SI_UNITX_SETTINGS = raw"
  \sisetup{
    detect-all,
    input-symbols = {()},
    table-format=2,
    table-number-alignment = center,
    input-open-uncertainty  = ,
    input-close-uncertainty = ,            
    table-figures-integer = 1,
    table-figures-decimal = 2,
    table-space-text-post = {\superscript{*}},
    table-space-text-post = {\superscript{**}},
    table-space-text-post = {\superscript{***}},
    table-space-text-pre    = (,
    table-space-text-post   = ),
  } % centering in tables
"


#textable
#=    Creates a well-formed Latex table from string components
  The content matrix consists of an array of matrices, which are overlayed
  via a one line offset. For example, the first matrix might be the values,
  the second matrix the errors. Observational rows can be put in
  desccontent, which is a vector of rows. Similarly, colnames is a vector
  of rows. So that column headers can span multiple columns, the widthcolnames
  and widthdesccontent provides the option to specify dimentions
  IN: colnames, rownames, content, desccontent (rows of
      descriptive data), notes, optionally widthcolnames and widthDescConent
      (which are the width of each entry in columns)
  OUT: The tex table as a string=#

function textable(;
  colnames::Vector{Vector{String}} = Vector{Vector{String}}(),
  contentrownames::Vector{String} = Vector{String}(),
  content::Vector{Matrix{String}} = Vector{Matrix{String}}(),
  descrownames::Vector{String} = Vector{String}(),
  desccontent::Vector{Vector{String}}=Vector{Vector{String}}(),
  notes::Vector{String} = Vector{String}(),
  linespacer::String = "\\\\",
  summarymathmode::Bool = true,
  widthcolnames::Vector{Vector{Int}} = #contains the number of columns for each entry
    broadcast((i::Int)->ones(Int,length(colnames[i])),1:length(colnames)),
  alignmentcolnames::Vector{Vector{String}} = #contains the number of columns for each entry
    broadcast((i::Int)->["r" for i ∈ 1:length(colnames[i])],1:length(colnames)),
  widthdesccontent::Vector{Vector{Int}} = #contains the number of columns for each entry
    broadcast((i::Int)->ones(Int,length(desccontent[i])),1:length(desccontent)),
  colheadername::Vector{String} = ["" for i ∈ 1:length(colnames)],
  alignmentstring = string(" l | ", join(["S" for i ∈ 1:length(desccontent[1])])),
  rowlabelheader::Bool = false,
  siunitxpreamble=true,
  siunitxsettings=SI_UNITX_SETTINGS)

  #size parameters for content
  numcontentrows::Int = length(contentrownames)

  numcontentcols::Int = sum(widthcolnames[1]) - rowlabelheader

  numcontentsubrows::Int = length(content)

  #intiate the stream
  b::IOBuffer = IOBuffer()

  siunitxpreamble && write(b,siunitxsettings)
  write(b, """
      \\begin{tabular}{$alignmentstring""")

    write(b,"}\n \\toprule")#filler tex

  for r::Int ∈ 1:length(colnames) #for each row of column headings
      #write(b,"\n\\\\[-1.8ex]") #use this if the below doesn't work
    write(b, colheadername[r])
    if length(colnames[r]) ≠ numcontentcols
      for c::Int ∈ 1:length(colnames[r]) #for each column heading
        (c≠1 || (!rowlabelheader)) && write(b,"\t&")
        write(b,"\t\\multicolumn{$(widthcolnames[r][c])}{$(alignmentcolnames[r][c])}{$(colnames[r][c])}")
      end
    else
      for c::Int ∈ 1:length(colnames[r]) #for each column heading
        (c≠1 || (!rowlabelheader)) && write(b,"\t&")
        write(b,"\t$(colnames[r][c])")
      end
    end
      write(b,"\n\\\\")
  end
  write(b, " \\midrule\n ")

  #now write out the table content
    if numcontentrows > 0
        write(b," \t \t ")
        for r::Int ∈ 1:numcontentrows
          write(b, "$(contentrownames[r])") #the row label
          for s::Int ∈ 1:numcontentsubrows #print the sub-rows for each row
            for c::Int ∈ 1:numcontentcols #print the columns for each sub-row
              write(b,"\t&\t$((content[s])[r, c])")
            end
            write(b,"\n $(linespacer) \t\t") #line-break stylistic formatting
          end
        end
        write(b, "\n \\midrule\n")
    end
    leftsummathflag::String = summarymathmode ? "" : "{"
    rightsummathflag::String = summarymathmode ? "" : "}"

    if length(desccontent)>0
        for r::Int ∈ 1:length(desccontent) #for each description row
          write(b, "$(descrownames[r])") #the row label
          if length(widthdesccontent[r]) ≠ numcontentcols
            for c::Int ∈ 1:length(desccontent[r]) #for each descriptive row
              if length(desccontent[r][c]) ≥ 1
                write(b,"\t&\t\\multicolumn{$(widthdesccontent[r][c])}{r}{$(leftsummathflag)$(desccontent[r][c])$(rightsummathflag)}")
              else
                write(b,"\t&\t\\multicolumn{$(widthdesccontent[r][c])}{r}{}")
              end
            end
          else
            for c::Int ∈ 1:length(desccontent[r]) #for each descriptive row
              if length(desccontent[r][c]) ≥ 1
                write(b,"\t&\t$(leftsummathflag)$(desccontent[r][c])$(rightsummathflag)")
              else
                write(b,"\t&\t")
              end
            end
          end
            write(b,"\n $(linespacer) \t\t") #line-break stylistic formatting
        end
    end

      write(b,""" \\bottomrule""")
  #if we have footnotes
  if length(notes) > 0
      write(b,"""\n \\\\[-1.0ex] \\textit{Notes:} \t """)
      for r ∈ length(notes)
          write(b, "\t&\t \\multicolumn{$numcontentcols}{l}{$(notes[r])}\n \\\\")
      end
  end
  write(b, "\t \\end{tabular}\n")
  return String(take!(copy(b)))
end


#=this version formats a table from a set of linear models
While customizable, the idea here is to have sensible default arguments
where possible.
IN: A collection of models, the rows to capture (returns an empty string in
the table if a row is not present, optionally column names, their width,
names of the rows, names of the descriptive rows, descriptive content,
a boolean switch for including signficance stars, signficance levels for the stars,
a scaling factor, the number of digits in values, and customizable
notes
OUT: A latex table string
=#
function textable(models::Vector{<:FMLM},
    getΣ::Union{Function, Vector{Function}},
    rows::Vector{String};
    colnames::Vector{Vector{String}} = [["{($i)}" for i ∈ 1:length(models)]],
    contentrownames::Vector{String} = String.(rows),
    descrownames::Vector{String}=Vector{String}(),
    desccontent::Vector{Vector{String}}=Vector{Vector{String}}(),
    notes::Vector{String} = Vector{String}(),
    linespacer::String = "\\\\",
    summarymathmode::Bool = true,
    widthcolnames::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(colnames[i])),1:length(colnames)),
    alignmentcolnames::Vector{Vector{String}} = #contains the number of columns for each entry
      broadcast((i::Int)->["r" for i ∈ 1:length(colnames[i])],1:length(colnames)),
    widthdesccontent::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(desccontent[i])),1:length(desccontent)),
    stars::Bool=true,
    starlvls::Vector{Float64} = [.9, .95, .99],
    starlegend::String = stars ? DEFAULT_STAR_LEGEND : "",
    starstrings =[
      raw"\textnormal{\superscript{*}}",
      raw"\textnormal{\superscript{**}}",
      raw"\textnormal{\superscript{***}}",],
    scaling::Vector{<:Real}=ones(length(rows)),
    decimaldigits::Int = 2,
    colheadername::Vector{String} = ["" for i::Int ∈ 1:length(colnames)],
    #alignmentstring::String = string(" l | ", join(["r" for i ∈ 1:length(desccontent[1])])),
    alignmentstring = string(" l | ", join(["S" for i ∈ 1:length(desccontent[1])])),
    rowlabelheader::Bool = false,
    pmethod=:ttest_2tailed,
    kwargs...)

  Ncols = length(models)
  numcontentrows = length(rows)

  #Pre-allocate the vector of SE errors
  σs::Vector{Vector{Float64}} =
      [Vector{Float64}(undef, models[i].K) for i∈1:Ncols]

  #pull out the β coefficients and N
  βs::Vector{Vector{Float64}} = [models[i].β for i::Int ∈ 1:Ncols]
  Ns::Vector{Int} = [models[i].N for i ∈ 1:Ncols]
  Xnames::Vector{Vector{String}} = [models[i].Xnames for i::Int ∈ 1:Ncols]

  #get the standard errors
  for c ∈ 1:Ncols
      if typeof(getΣ) <: Function
        σs[c] .= sqrt.(diag(getΣ(models[c])))
      else
        σs[c] .= sqrt.(diag(getΣ[c](models[c])))
      end
  end

  if length(starlegend) > 0
      notes = [starlegend; notes]
  end

  #get the content matrices
  content::Vector{Matrix{String}} = getcontentmatrices!( βs, σs,
    Xnames, Ns, rows,
    stars, starlvls, scaling,
    decimaldigits, starstrings,pmethod)

  return textable(;colnames,
    contentrownames, content, descrownames,
    desccontent, notes, linespacer,
    summarymathmode, widthcolnames,
    widthdesccontent, colheadername,
    alignmentcolnames,
    alignmentstring, rowlabelheader, kwargs...)
end



calculatep(::Val{:ttest_2tailed};N,β,σ) = cdf(TDist(N), β/σ) > .5 ? 1-(1.0 - cdf(TDist(N), β/σ))*2.0 : 1.0 - cdf(TDist(N), β/σ)*2.0
calculatep(::Val{:ttest_2tailed_offset};N,β,σ,β0) = calculatep(Val{:ttest_2tailed}();N,β-β0,σ)
calculatep(::Val{:ttest_1tailed};N,β,σ) = cdf(TDist(N), β/σ) > .5 ? 1-(1.0 - cdf(TDist(N), β/σ)) : 1.0 - cdf(TDist(N), β/σ)
calculatep(::Val{:asymptotic_2tailed};N,β,σ) = cdf(Normal(), β/σ) > .5 ? 1-(1.0 - cdf(Normal(), β/σ))*2.0 : 1.0 - cdf(Normal(), β/σ)*2.0
calculatep(::Val{:asymptotic_1tailed};N,β,σ) = cdf(Normal(), β/σ) > .5 ? 1-(1.0 - cdf(Normal(), β/σ)) : 1.0 - cdf(Normal(), β/σ)


#=getcontentmatrices!
this is a bit of utiltiy code for making tex tables. It is not model specific
hence why it was extracted into its own function
#IN: a vector of βs, σs, names, rows (coefficeints selected). Optional parameters
include a pre-allcoation of the string matrix, a switch for the inclusion of
stars, a scaling factor, and the number of digits (rounding level)=#
#OUT: Writes to and returns the content matrix
function getcontentmatrices!(;
    βs::Vector{Vector{Float64}},
    σs::Vector{Vector{Float64}},
    Xnames::Vector{Vector{String}},
    Ns::Vector{Int},
    rows::Vector{String},
    psasσs::Bool = false,

    stars::Bool=true, #whether to display signficance stars
    starlvls::Vector{Float64} = [.9, .95, .99],  #cutoffs for signficance (must be sorted)
    #starstrings::Vector{String} =
    #  ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    starstrings =[
      raw"\textnormal{\superscript{*}}",
      raw"\textnormal{\superscript{**}}",
      raw"\textnormal{\superscript{***}}",],
    scaling::Vector{<:Real}=ones(length(rows)), # an optional scaling factor
    decimaldigits::Int = 2,
    pmethod=:ttest_2tailed) #number of decimal digits

    K= length(βs)

    #the below allows for customization of the p value star techniques
    #can either supply a singleton or a vector for all
    getpfunc(m::Val) = (;kwargs...)->calculatep(m; kwargs...)
    getpfunc(f::Function) = f
    getpfunc(v::Symbol) = getpfunc(Val(v))
    getpfunc(v::Vector) = v .|> getpfunc
    pfuncs = getpfunc(pmethod)

    content::Vector{Matrix{String}} = #will hold the coefficients and errors
        [fill("",length(rows),length(βs)) for i ∈ 1:2];

    #iterate through the models
    for c ∈ 1:K

        #build a dictionary of the names
        xnametable::Dict = Dict(Xnames[c][i] => i for i ∈ 1:length(βs[c]))

        for r ∈ 1:length(rows)
            if haskey(xnametable, rows[r]) #need to check if it exists

              if stars
                ind = xnametable[rows[r]]

                if psasσs
                  p = σs[c][ind]
                else
                  p = pfuncs[c](N=Ns[c], β=βs[c][ind], σ=σs[c][ind])
                  #p = cdf(TDist(Ns[c]), βs[c][ind]/σs[c][ind]) #get CDF from T distribution
                  #p = p > .5 ? 1-(1.0 - p)*2.0 : 1.0 - p*2.0 #calc 2-tailed p value
                end
                starstring = sum(p.>starlvls) > 0 ?  starstrings[sum(p.>starlvls)] : ""
              else
                starstring = ""
              end

              #scale, round and write the β coefficeint and σ into the string matrices
              content[1][r,c] =
                #"\$$(num2str(βs[c][ind], decimaldigits, scalefactor=scaling[r]))\$$starstring"
                "$(num2str(βs[c][ind], decimaldigits, scalefactor=scaling[r]))$starstring"
              content[2][r,c] =
                #"(\$$(num2str(σs[c][ind], decimaldigits, scalefactor=scaling[r]))\$)"
                "($(num2str(σs[c][ind], decimaldigits, scalefactor=scaling[r])))"
          end
        end
    end

  return content
end

#=Writes the fully formed string to the file
  IN: A string with the output, optionally a path to the file, and the
  output name
  OUT: Writes the string to a file=#
function writetextable(tablestring::String;
  path::String=pwd(),
  outname::String = "$path/results$(Dates.format(now(),"yymmdd_THMS")).tex")

  oStream::IOStream = open("$path/$outname","w+")
  write(oStream, tablestring)
  close(oStream)
end


function iotest()
  #test parameters
  nrowsContent = 8
  ncols = 4
  nsecs = ncols ÷ 2
  nrowsDesc = 2
  nsubRows = 2

  path = pwd() * "/results"
  outname = "test.tex"
  footerName = "footer.tex"
  headerName = "header.tex"

  #test table
  colnames::Vector{Vector{String}} =
      [broadcast(i->"secs$i",1:nsecs),broadcast(i->"ncols$i",1:ncols)]
  contentrownames::Vector{String} = ["nrows$i" for i ∈ 1:nrowsContent]
  content::Vector{Matrix{String}} =
      (string).(broadcast((i::Int)->i .* ones(nrowsContent,ncols), 1:nsubRows))
  descrownames::Vector{String} = ["Desc. row $i" for i ∈ 1:nrowsDesc]
  desccontent::Vector{Vector{String}} =
      [broadcast(i->"desc-r1c$i",1:ncols),broadcast(i->"desc-r1c$i",1:ncols)]
  notes::Vector{String} = ["a note"]
  widthcolnames::Vector{Vector{Int}} =
      [broadcast(i->(ncols ÷ nsecs), 1:nsecs),broadcast(i->1, 1:ncols)]

  s=textable(;colnames,  contentrownames, content, descrownames,
      desccontent, notes, widthcolnames)

  println(array2string(fill(.45/π,3,3,4),decimaldigits=4))

  print(s)
end
