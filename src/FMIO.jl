
######################IO################


#functions for converting Floats to Strings
function num2Str(x::T, decimals::Int=DEFAULT_DECIMALS;
    scaleFactor::Float64 = 1.0,
    scaleFactorOnHurdle::Float64 = 1.0,
    Ints::Bool = false, scaleHurdle::U=nothing, texEquation::Bool = false)::String where {T<:Real, U<:Union{Nothing,Real}}

  local outStr::String

  x*=scaleFactor

  if scaleHurdle ≠ nothing && x≥scaleHurdle
    outStr = "$(Int(round(x*scaleFactorOnHurdle)))"
  elseif decimals == 0 || (Ints && round(x) == x)
    outStr = "$(Int(round(x)))"
  else
    outStr = format("{:.$(decimals)f}",x)
  end

  outStr = texEquation ? "\$$outStr\$" : outStr

  return outStr
end

#converts a vector of anything to a string given an abstract delimitor
function vec2String(v::Vector{T}, delim::W = "")::String where
  {T<:Any, W<:Union{AbstractString,Char}}

  return join(String.(v), delim)
end



#texTable
#=    Creates a well-formed Latex table from string components
  The content matrix consists of an array of matrices, which are overlayed
  via a one line offset. For example, the first matrix might be the values,
  the second matrix the errors. Observational rows can be put in
  descContent, which is a vector of rows. Similarly, colNames is a vector
  of rows. So that column headers can span multiple columns, the widthColNames
  and widthDescContent provides the option to specify dimentions
  IN: caption (title), colNames, rowNames, content, descContent (rows of
      descriptive data), notes, optionally widthColNames and widthDescConent
      (which are the width of each entry in columns)
  OUT: The tex table as a string=#

function texTable(; titleCaption::String = "",
  caption::String = "",
  colNames::Vector{Vector{String}} = Vector{Vector{String}}(),
  contentRowNames::Vector{String} = Vector{String}(),
  content::Vector{Matrix{String}} = Vector{Matrix{String}}(),
  descRowNames::Vector{String} = Vector{String}(),
  descContent::Vector{Vector{String}}=Vector{Vector{String}}(),
  notes::Vector{String} = Vector{String}(),
  arrayStretch::Float64 = 1.5,
  lineSpacer::String = "\\\\",
  summaryMathMode::Bool = true,
  widthColNames::Vector{Vector{Int}} = #contains the number of columns for each entry
    broadcast((i::Int)->ones(Int,length(colNames[i])),1:length(colNames)),
  alignmentColNames::Vector{Vector{String}} = #contains the number of columns for each entry
    broadcast((i::Int)->["r" for i ∈ 1:length(colNames[i])],1:length(colNames)),
  widthDescContent::Vector{Vector{Int}} = #contains the number of columns for each entry
    broadcast((i::Int)->ones(Int,length(descContent[i])),1:length(descContent)),
  colHeaderName::Vector{String} = ["" for i ∈ 1:length(colNames)],
  nakedTable::Bool = false)

  #size parameters for content
  numContentRows::Int = length(contentRowNames)

  numContentCols::Int = sum(widthColNames[1])

  numContentSubRows::Int = length(content)

      #intiate the stream
  b::IOBuffer = IOBuffer()
  (!nakedTable) && write(b, """
      %This table was programatically generated

      \\begin{table} \\caption{$titleCaption} \\label{} \\centering
      $(length(caption)>0 ? "\\textit{$caption\\\\}" : "")
      \\renewcommand{\\arraystretch}{$arrayStretch}""")
  write(b, """
      \\begin{tabular}{l | """)

    for i ∈ 1:(numContentCols) #set the dimensions in tex
      write(b,"r")
    end
    #write(b, "{\\textwidth}{Xccccccc}")
    write(b,"}\n \\toprule")#filler tex

  for r::Int ∈ 1:length(colNames) #for each row of column headings
      #write(b,"\n\\\\[-1.8ex]") #use this if the below doesn't work
    write(b, colHeaderName[r])
    if length(colNames[r]) ≠ numContentCols
      for c::Int ∈ 1:length(colNames[r]) #for each column heading
          write(b,"\t&\t\\multicolumn{$(widthColNames[r][c])}{$(alignmentColNames[r][c])}{$(colNames[r][c])}")
      end
    else
      for c::Int ∈ 1:length(colNames[r]) #for each column heading
        write(b,"\t&\t$(colNames[r][c])")
      end
    end
      write(b,"\n\\\\")
  end
  write(b, " \\midrule\n ")

  #now write out the table content
    if numContentRows > 0
        write(b," \t \t ")
        for r::Int ∈ 1:numContentRows
          write(b, "$(contentRowNames[r])") #the row label
          for s::Int ∈ 1:numContentSubRows #print the sub-rows for each row
            for c::Int ∈ 1:numContentCols #print the columns for each sub-row
              #write(b,"\t&\t\\multicolumn{1}{r}{\$$((content[s])[r, c])\$}")
              write(b,"\t&\t$((content[s])[r, c])")
            end
            write(b,"\n $(lineSpacer) \t\t") #line-break stylistic formatting
          end
        end
        write(b, "\n \\midrule\n")
    end
    sumMathFlag::String = summaryMathMode ? "\$" : ""

    if length(descContent)>0
        for r::Int ∈ 1:length(descContent) #for each description row
          write(b, "$(descRowNames[r])") #the row label
          if length(widthDescContent[r]) ≠ numContentCols
            for c::Int ∈ 1:length(descContent[r]) #for each descriptive row
              if length(descContent[r][c]) ≥ 1
                write(b,"\t&\t\\multicolumn{$(widthDescContent[r][c])}{r}{$(sumMathFlag)$(descContent[r][c])$(sumMathFlag)}")
              else
                write(b,"\t&\t\\multicolumn{$(widthDescContent[r][c])}{r}{}")
              end
            end
          else
            for c::Int ∈ 1:length(descContent[r]) #for each descriptive row
              if length(descContent[r][c]) ≥ 1
                write(b,"\t&\t$(sumMathFlag)$(descContent[r][c])$(sumMathFlag)")
              else
                write(b,"\t&\t")
              end
            end
          end
            write(b,"\n $(lineSpacer) \t\t") #line-break stylistic formatting
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
  (!nakedTable) && write(b, " \\\\ \\end{table}")
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
function texTable(models::Vector{FMLM}, getΣ::Function, rows::Vector{Symbol};
    titleCaption::String = "",
    caption::String = "",
    colNames::Vector{Vector{String}} = [["C$i" for i ∈ 1:length(models)]],
    contentRowNames::Vector{String} = String.(rows),
    descRowNames::Vector{String}=Vector{String}(),
    descContent::Vector{Vector{String}}=Vector{Vector{String}}(),
    notes::Vector{String} = Vector{String}(),
    arrayStretch::Float64 = 1.5,
    lineSpacer::String = "\\\\",
    summaryMathMode::Bool = true,
    widthColNames::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(colNames[i])),1:length(colNames)),
    alignmentColNames::Vector{Vector{String}} = #contains the number of columns for each entry
      broadcast((i::Int)->["r" for i ∈ 1:length(colNames[i])],1:length(colNames)),
    widthDescContent::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(descContent[i])),1:length(descContent)),
    stars::Bool=true,
    starLvls::Vector{Float64} = [.9, .95, .99],
    starLegend::String = stars ? DEFAULT_STAR_LEGEND : "",
    starStrings::Vector{String} =
      ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    scaling::Vector{Float64}=ones(length(rows)),
    decimalDigits::Int = 2,
    colHeaderName::Vector{String} = ["" for i::Int ∈ 1:length(colNames)],
    nakedTable::Bool = false)

  numCols = length(models)
  numContentRows = length(rows)

  #Pre-allocate the vector of SE errors
  modelsσ::Vector{Vector{Float64}} =
      [Vector{Float64}(undef, models[i].K) for i∈1:numCols]

  #pull out the β coefficients and N
  modelsβ::Vector{Vector{Float64}} = [models[i].β for i::Int ∈ 1:numCols]
  modelsN::Vector{Int} = [models[i].N for i ∈ 1:numCols]
  modelsXNames::Vector{Vector{Symbol}} = [models[i].XNames for i::Int ∈ 1:numCols]

  #get the standard errors
  for c ∈ 1:numCols
      modelsσ[c] .= sqrt.(diag(getΣ(models[c])))
  end

  if length(starLegend) > 0
      notes = [starLegend; notes]
  end

  #get the content matrices
  content::Vector{Matrix{String}} = getContentMatrices!( modelsβ, modelsσ, modelsXNames, modelsN, rows,
      stars=stars, starLvls=starLvls, scaling=scaling,
      decimalDigits=decimalDigits, starStrings=starStrings)

  return texTable(titleCaption = titleCaption,caption=caption, colNames=colNames,
    contentRowNames=contentRowNames, content=content, descRowNames=descRowNames,
    descContent=descContent, notes=notes, arrayStretch=arrayStretch, lineSpacer=lineSpacer,
    summaryMathMode=summaryMathMode, widthColNames=widthColNames,
    widthDescContent=widthDescContent, colHeaderName=colHeaderName,
    alignmentColNames=alignmentColNames, nakedTable=nakedTable)
end

function texTable(models::Vector{FM2SLS}, getΣ::Function, rows::Vector{Symbol};
    titleCaption::String = "",
    caption::String = "",
    colNames::Vector{Vector{String}} = [["C$i" for i ∈ 1:length(models)]],
    contentRowNames::Vector{String} = String.(rows),
    descRowNames::Vector{String}=Vector{String}(),
    descContent::Vector{Vector{String}}=Vector{Vector{String}}(),
    notes::Vector{String} = Vector{String}(),
    arrayStretch::Float64 = 1.5,
    lineSpacer::String = "\\\\",
    summaryMathMode::Bool = true,
    widthColNames::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(colNames[i])),1:length(colNames)),
    alignmentColNames::Vector{Vector{String}} = #contains the number of columns for each entry
      broadcast((i::Int)->["r" for i ∈ 1:length(colNames[i])],1:length(colNames)),
    widthDescContent::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(descContent[i])),1:length(descContent)),
    stars::Bool=true,
    starLvls::Vector{Float64} = [.9, .95, .99],
    starLegend::String = stars ? DEFAULT_STAR_LEGEND : "",
    starStrings::Vector{String} =
      ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    scaling::Vector{Float64}=ones(length(rows)),
    decimalDigits::Int = 2,
    colHeaderName::Vector{String} = ["" for i ∈ 1:length(colNames)],
    nakedTable=nakedTable)

    numCols = length(models)
    numContentRows = length(rows)

    #Pre-allocate the vector of SE errors
    modelsσ::Vector{Vector{Float64}} =
    [Vector{Float64}(models[i].K+models[i].KW) for i∈1:numCols]


    #pull out the β coefficients and N
    modelsβ::Vector{Vector{Float64}} = [models[i].δ2 for i ∈ 1:numCols]
    modelsN::Vector{Int} = [models[i].N for i ∈ 1:numCols]
    modelsXWNames::Vector{Vector{Symbol}} = [[models[i].XNames; models[i].WNames]  for i ∈ 1:numCols]

    #get the standard errors
    for c ∈ 1:numCols
        modelsσ[c] .= sqrt.(diag(getΣ(models[c])))
    end

    if length(starLegend) > 0
        notes = [starLegend; notes]
    end

    #get the content matrices
    content::Vector{Matrix{String}} = getContentMatrices!( modelsβ, modelsσ, modelsXWNames, modelsN, rows,
      stars=stars, starLvls=starLvls, scaling=scaling,
      decimalDigits=decimalDigits, starStrings=starStrings)

    texTable(titleCaption = titleCaption,caption=caption, colNames=colNames,
      contentRowNames=contentRowNames, content=content, descRowNames=descRowNames,
      descContent=descContent, notes=notes, arrayStretch=arrayStretch, lineSpacer=lineSpacer,
      summaryMathMode=summaryMathMode, widthColNames=widthColNames,
      widthDescContent=widthDescContent, colHeaderName=colHeaderName,
      alignmentColNames=alignmentColNames, nakedTable=nakedTable)
end



#=getContentMatrices!
this is a bit of utiltiy code for making tex tables. It is not model specific
hence why it was extracted into its own function
#IN: a vector of βs, σs, names, rows (coefficeints selected). Optional parameters
include a pre-allcoation of the string matrix, a switch for the inclusion of
stars, a scaling factor, and the number of digits (rounding level)=#
#OUT: Writes to and returns the content matrix
function getContentMatrices!(modelsβ::Vector{Vector{Float64}},
    modelsσ::Vector{Vector{Float64}},
    modelsXNames::Vector{Vector{Symbol}},
    modelsN::Vector{Int},
    rows::Vector{Symbol};
    stars::Bool=true, #whether to display signficance stars
    starLvls::Vector{Float64} = [.9, .95, .99],  #cutoffs for signficance (must be sorted)
    starStrings::Vector{String} =
      ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    scaling::Vector{Float64}=ones(length(rows)), # an optional scaling factor
    decimalDigits::Int = 2) #number of decimal digits

    content::Vector{Matrix{String}} = #will hold the coefficients and errors
        [fill("",length(rows),length(modelsβ)) for i ∈ 1:2];

    #iterate through the models
    for c ∈ 1:length(modelsβ)

        #build a dictionary of the names
        XNameTbl::Dict = Dict(modelsXNames[c][i] => i for i::Int ∈ 1:length(modelsβ[c]))

        for r ∈ 1:length(rows)
            if haskey(XNameTbl, rows[r]) #need to check if it exists
                ind::Int = XNameTbl[rows[r]]
                p::Float64 =
                    cdf(TDist(modelsN[c]), modelsβ[c][ind]/modelsσ[c][ind]) #get CDF from T distribution
                p = p > .5 ? 1-(1.0 - p)*2.0 : 1.0 - p*2.0 #calc 2-tailed p value
                sigLevel::Int = sum(p.>starLvls)
                if sigLevel > 0 && stars
                    starString::String = starStrings[sum(p.>starLvls)]
                else
                    starString = ""
                end

                #scale, round and write the β coefficeint and σ into the string matrices
                content[1][r,c] =
                  #"$(round(scaling[r]*modelsβ[c][ind],decimalDigits))^{$starString}"
                  "\$$(num2Str(modelsβ[c][ind], decimalDigits, scaleFactor=scaling[r]))\$$starString"
                content[2][r,c] =
                  "(\$$(num2Str(modelsσ[c][ind], decimalDigits, scaleFactor=scaling[r]))\$)"
                  #"($(round(scaling[r]*modelsσ[c][ind],decimalDigits)))"
            end
        end
    end

  return content
end

#=Writes the fully formed string to the file
  IN: A string with the output, optionally a path to the file, and the
  output name
  OUT: Writes the string to a file=#
function writeTables2File(tablesString::String;
    path::String=pwd(),
    outName::String = "$path\\results$(Dates.format(now(),"yymmdd_THMS")).tex")

    oStream::IOStream = open("$path\\$outName","w+")
    write(oStream, tablesString)
    close(oStream)
end

#=Writes the table strings to a file
  Convenience method which takes an array of table strings and optionally
  header and footer strings and writes the table to a tex file
  IN: A string array of table strings, optionally a header string, a footer
  string, a path to the file, the output name
  OUT: Writes the tables to a file=#
function writeTables2File(tables::Vector{String};
    header::String = "\n", footer::String = "",
    path::String=pwd(),
    outName::String = "$path\\results$(Dates.format(now(),"yymmdd_THMS")).tex")

    b::IOBuffer = IOBuffer() #make an IO buffer
    write(b,"\n$header\n\n")
    for t ∈ tables
        write(b,"$t\n\n\n") #write each table
    end
    write(b,footer) #write the footer

    #write the output via another instance
    writeTables2File(String(take!(b)), path=path, outName = outName)
end

#=Writes the table strings to a file
Convenience method which takes an array of table strings and the locations of
header and footer files.
IN: A string array of table strings, the name of a header file, the name of a
a footer file, optionally a path to the file, the output name
OUT: Writes the tables to a file=#
function writeTables2File(tables::Vector{String}, headerName::String,
  footerName::String;
  path::String=pwd(),
  outName::String = "$path\\results$(Dates.format(now(),"yymmdd_THMS")).tex")

  iStream::IOStream = open("$path\\$headerName")
  header::String = String(take!(IOBuffer(read(iStream))))
  close(iStream)

  iStream = open("$path\\$footerName")
  footer::String = String(take!(IOBuffer(read(iStream))))
  close(iStream)

  #write the output via another instance
  writeTables2File(tables, path=path, outName = outName, header=header,
      footer=footer)
end

function writeNakedTable(tableString::String;
  path::String=pwd(),
  outName::String = "$path\\results$(Dates.format(now(),"yymmdd_THMS")).tex")

  oStream::IOStream = open("$path\\$outName","w+")
  write(oStream, tableString)
  close(oStream)
end

#array2String
  #=takes content of an array, formats each cell into a string, appends a
  prefix and suffix, and sends back an array of strings with the same
  dimension as the input
  IN: A general array, optionally a prefix, suffix, scaling factor, and
   set number of decimal digits
  OUT: An array of strings with the same dimensions as the input=#
function array2String(m::Array{Float64}; prefix::String = "",
  suffix::String = "", scaling::Float64 = 1.0,
  decimalDigits::Int = 2)::Array{String}

  #broadcast the prefix, value in m, and suffix and concatenate respectively
  stringMat::Array{String} =
      ((x::Float64)->(prefix*"$(round(x,decimalDigits))"*suffix)).(m)

  return stringMat
end

#vcats like with rbind in r
function vbind(df1::AbstractDataFrame, df2::AbstractDataFrame)::AbstractDataFrame
  #get info on what we are joining
  fields1::Vector{Symbol} = names(df1)
  fields2::Vector{Symbol} = names(df2)
  allFields::Vector{Symbol} = union(fields1, fields2)

  rows1::Int = size(df1,1)
  rows2::Int = size(df2,1)

  for f::Symbol ∈ setdiff(allFields, fields1)
    origType::Type = eltype(df2[f])

    if !(Missing <: origType) #in this case we need to promote the eltype
      origType = Union{origType, Missing}
      df2[f] = Vector{origType}(df2[f])
    end

    df1[ f] = Vector{origType}(missing, rows1)
  end

  for f::Symbol ∈ setdiff(allFields, fields2)
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




function IOTest()
  #test parameters
  nrowsContent = 8
  ncols = 4
  nsecs = ncols ÷ 2
  nrowsDesc = 2
  nsubRows = 2

  path = pwd() * "\\results"
  outName = "test.tex"
  footerName = "footer.tex"
  headerName = "header.tex"

  #test table
  caption = "test table"
  colNames::Vector{Vector{String}} =
      [broadcast(i->"secs$i",1:nsecs),broadcast(i->"ncols$i",1:ncols)]
  contentRowNames::Vector{String} = ["nrows$i" for i ∈ 1:nrowsContent]
  content::Vector{Matrix{String}} =
      array2String.(broadcast((i::Int)->i .* ones(nrowsContent,ncols), 1:nsubRows))
  descRowNames::Vector{String} = ["Desc. row $i" for i ∈ 1:nrowsDesc]
  descContent::Vector{Vector{String}} =
      [broadcast(i->"desc-r1c$i",1:ncols),broadcast(i->"desc-r1c$i",1:ncols)]
  notes::Vector{String} = ["a note"]
  widthColNames::Vector{Vector{Int}} =
      [broadcast(i->(ncols ÷ nsecs), 1:nsecs),broadcast(i->1, 1:ncols)]

  s=texTable(caption, colNames,  contentRowNames, content, descRowNames,
      descContent, notes, widthColNames=widthColNames)

  writeTables2File([s,s], headerName, footerName, path=path, outName = outName)

  #println(array2String(fill(.45/π,3,3,4),decimalDigits=4))

  print(s)
end
