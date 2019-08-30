######################IO################

vec2String(args...; keyargs...) = error("Depreciated- replace with join")
texTable(args...; keyargs...) = error("Depreciated- replace with textable")
writeTexTable(args...; keyargs...) = error("Depreciated- replace with writetextable")
array2String(args...; keyargs...) = error("Depreciated- replace with array2string")

#functions for converting Floats to Strings
function num2str(x::T, decimals::Int=DEFAULT_DECIMALS;
    scalefactor::Float64 = 1.0,
    scalefactoronhurdle::Float64 = 1.0,
    Ints::Bool = false, scalehurdle::U=nothing, texequation::Bool = false)::String where {
      T<:Union{Real,Missing,Nothing}, U<:Union{Nothing,Real}}

  local outStr::String

  x*=scalefactor

  #this controls for both missing and nothing values
  if ismissing(something(x,missing))
    outStr = ""
  elseif scalehurdle ≠ nothing && x≥scalehurdle
    outStr = "$(Int(round(x*scalefactoronhurdle)))"
  elseif decimals == 0 || (Ints && round(x) == x)
    outStr = "$(Int(round(x)))"
  else
    outStr = format("{:.$(decimals)f}",x)
  end

  outStr = texequation ? "\$$outStr\$" : outStr

  return outStr
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
  alignmentstring::String = string(" l | ", join(["r" for i ∈ 1:length(desccontent[1])])))

  #size parameters for content
  numContentRows::Int = length(contentrownames)

  numContentCols::Int = sum(widthcolnames[1])

  numContentSubRows::Int = length(content)

      #intiate the stream
  b::IOBuffer = IOBuffer()
  write(b, """
      \\begin{tabular}{$alignmentstring""")

    write(b,"}\n \\toprule")#filler tex

  for r::Int ∈ 1:length(colnames) #for each row of column headings
      #write(b,"\n\\\\[-1.8ex]") #use this if the below doesn't work
    write(b, colheadername[r])
    if length(colnames[r]) ≠ numContentCols
      for c::Int ∈ 1:length(colnames[r]) #for each column heading
          write(b,"\t&\t\\multicolumn{$(widthcolnames[r][c])}{$(alignmentcolnames[r][c])}{$(colnames[r][c])}")
      end
    else
      for c::Int ∈ 1:length(colnames[r]) #for each column heading
        write(b,"\t&\t$(colnames[r][c])")
      end
    end
      write(b,"\n\\\\")
  end
  write(b, " \\midrule\n ")

  #now write out the table content
    if numContentRows > 0
        write(b," \t \t ")
        for r::Int ∈ 1:numContentRows
          write(b, "$(contentrownames[r])") #the row label
          for s::Int ∈ 1:numContentSubRows #print the sub-rows for each row
            for c::Int ∈ 1:numContentCols #print the columns for each sub-row
              write(b,"\t&\t$((content[s])[r, c])")
            end
            write(b,"\n $(linespacer) \t\t") #line-break stylistic formatting
          end
        end
        write(b, "\n \\midrule\n")
    end
    sumMathFlag::String = summarymathmode ? "\$" : ""

    if length(desccontent)>0
        for r::Int ∈ 1:length(desccontent) #for each description row
          write(b, "$(descrownames[r])") #the row label
          if length(widthdesccontent[r]) ≠ numContentCols
            for c::Int ∈ 1:length(desccontent[r]) #for each descriptive row
              if length(desccontent[r][c]) ≥ 1
                write(b,"\t&\t\\multicolumn{$(widthdesccontent[r][c])}{r}{$(sumMathFlag)$(desccontent[r][c])$(sumMathFlag)}")
              else
                write(b,"\t&\t\\multicolumn{$(widthdesccontent[r][c])}{r}{}")
              end
            end
          else
            for c::Int ∈ 1:length(desccontent[r]) #for each descriptive row
              if length(desccontent[r][c]) ≥ 1
                write(b,"\t&\t$(sumMathFlag)$(desccontent[r][c])$(sumMathFlag)")
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
          write(b, "\t&\t \\multicolumn{$numContentCols}{l}{$(notes[r])}\n \\\\")
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
function textable(models::Vector{FMLM},
    getΣ::T where T<:Union{Function, Vector{Function}},
    rows::Vector{Symbol};
    colnames::Vector{Vector{String}} = [["($i)" for i ∈ 1:length(models)]],
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
    starstrings::Vector{String} =
      ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    scaling::Vector{Float64}=ones(length(rows)),
    decimaldigits::Int = 2,
    colheadername::Vector{String} = ["" for i::Int ∈ 1:length(colnames)],
    alignmentstring::String = string(" l | ", join(["r" for i ∈ 1:length(desccontent[1])])))

  Ncols = length(models)
  numContentRows = length(rows)

  #Pre-allocate the vector of SE errors
  σs::Vector{Vector{Float64}} =
      [Vector{Float64}(undef, models[i].K) for i∈1:Ncols]

  #pull out the β coefficients and N
  βs::Vector{Vector{Float64}} = [models[i].β for i::Int ∈ 1:Ncols]
  Ns::Vector{Int} = [models[i].N for i ∈ 1:Ncols]
  modelsxnames::Vector{Vector{Symbol}} = [models[i].xnames for i::Int ∈ 1:Ncols]

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
  content::Vector{Matrix{String}} = getcontentmatrices!( βs, σs, xnames, Ns, rows,
      stars=stars, starlvls=starlvls, scaling=scaling,
      decimaldigits=decimaldigits, starstrings=starstrings)

  return textable(colnames=colnames,
    contentrownames=contentrownames, content=content, descrownames=descrownames,
    desccontent=desccontent, notes=notes, linespacer=linespacer,
    summarymathmode=summarymathmode, widthcolnames=widthcolnames,
    widthdesccontent=widthdesccontent, colheadername=colheadername,
    alignmentcolnames=alignmentcolnames,
    alignmentstring=alignmentstring)
end


#=getcontentmatrices!
this is a bit of utiltiy code for making tex tables. It is not model specific
hence why it was extracted into its own function
#IN: a vector of βs, σs, names, rows (coefficeints selected). Optional parameters
include a pre-allcoation of the string matrix, a switch for the inclusion of
stars, a scaling factor, and the number of digits (rounding level)=#
#OUT: Writes to and returns the content matrix
function getcontentmatrices!(;
    βs::Vector{Vector{Float64}} = error("βs is a required argument for getcontentmatrices"),
    σs::Vector{Vector{Float64}} = error("σs is a required argument for getcontentmatrices"),
    xnames::Vector{Vector{Symbol}} = error("xnames is a required argument for getcontentmatrices"),
    Ns::Vector{Int} = error("Ns is a required argument for getcontentmatrices"),
    rows::Vector{Symbol} = error("rows is a required argument for getcontentmatrices"),
    stars::Bool=true, #whether to display signficance stars
    starlvls::Vector{Float64} = [.9, .95, .99],  #cutoffs for signficance (must be sorted)
    starstrings::Vector{String} =
      ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    scaling::Vector{Float64}=ones(length(rows)), # an optional scaling factor
    decimaldigits::Int = 2) #number of decimal digits

    content::Vector{Matrix{String}} = #will hold the coefficients and errors
        [fill("",length(rows),length(βs)) for i ∈ 1:2];

    #iterate through the models
    for c ∈ 1:length(βs)

        #build a dictionary of the names
        XNameTbl::Dict = Dict(xnames[c][i] => i for i::Int ∈ 1:length(βs[c]))

        for r ∈ 1:length(rows)
            if haskey(XNameTbl, rows[r]) #need to check if it exists
                ind::Int = XNameTbl[rows[r]]
                p::Float64 =
                    cdf(TDist(Ns[c]), βs[c][ind]/σs[c][ind]) #get CDF from T distribution
                p = p > .5 ? 1-(1.0 - p)*2.0 : 1.0 - p*2.0 #calc 2-tailed p value
                sigLevel::Int = sum(p.>starlvls)
                if sigLevel > 0 && stars
                    starString::String = starstrings[sum(p.>starlvls)]
                else
                    starString = ""
                end

                #scale, round and write the β coefficeint and σ into the string matrices
                content[1][r,c] =
                  #"$(round(scaling[r]*βs[c][ind],decimaldigits))^{$starString}"
                  "\$$(num2str(βs[c][ind], decimaldigits, scalefactor=scaling[r]))\$$starString"
                content[2][r,c] =
                  "(\$$(num2str(σs[c][ind], decimaldigits, scalefactor=scaling[r]))\$)"
                  #"($(round(scaling[r]*σs[c][ind],decimaldigits)))"
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
  outname::String = "$path\\results$(Dates.format(now(),"yymmdd_THMS")).tex")

  oStream::IOStream = open("$path\\$outname","w+")
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

  path = pwd() * "\\results"
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

  s=textable(colnames,  contentrownames, content, descrownames,
      desccontent, notes, widthcolnames=widthcolnames)

  println(array2string(fill(.45/π,3,3,4),decimaldigits=4))

  print(s)
end
