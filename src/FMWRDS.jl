#contains methods related to wrds
#for brevity of the finometrics file, long-string constants are kept here
using ODBC

#DSN_NAME = "wrds-pgdata-64"
WRDS_DSN = "Driver={PostgreSQL ANSI(x64)};Server=wrds-pgdata.wharton.upenn.edu;
         Port=9737;Database=wrds;Uid=ctepp111;Pwd=CSTWH1211@ce;sslmode=require;"

#= Default is
dsn = ODBC.DSN("Driver={PostgreSQL ANSI(x64)};Server=wrds-pgdata.wharton.upenn.edu;
         Port=9737;Database=wrds;Uid=ctepp111;Pwd=CSTWH1211@ce;sslmode=require;",
         prompt=false)
=#

QUICK_QUERIES = Dict(:libraries => "select distinct table_schema
                   from information_schema.tables
                   where table_type ='VIEW'
                   or table_type ='FOREIGN TABLE'
                   order by table_schema",
                   :crspTables => "select distinct table_name
                   from information_schema.columns
                   where table_schema='crsp'
                   order by table_name",
                   :crspaTables => "select distinct table_name
                   from information_schema.columns
                   where table_schema='crspa'
                   order by table_name",
                   :crspDailyMeta => "select column_name
                   from information_schema.columns
                   where table_schema='crsp'
                   and table_name='dsf'
                   order by column_name",
                   :crspMonthlyMeta => "select column_name
                   from information_schema.columns
                   where table_schema='crsp'
                   and table_name='msf'
                   order by column_name",
                   :exCrsp => "select cusip,permno,date,bidlo,askhi
                   from crsp.dsf
                   where permno in (14593, 90319, 12490, 17778)
                   and date between '2010-01-01' and '2013-12-31'
                   and askhi > 2000",
                   :sp500 => "select *
                   from crspa.msp500
                   where caldt between '1926-01-01' and '2018-12-31'"
                   )

function wrdsGet(connectString::String = WRDS_DSN)
  return ODBC.DSN(connectString, prompt=false)
end

function wrdsDrop!(wrds::ODBC.DSN)::Nothing
  ODBC.disconnect!(wrds)
  return nothing
end

function wrdsQuery(wrds::ODBC.DSN, qstring::String;
    limit::T where T<:Union{Nothing, Int} = nothing)::DataFrame

  local suffix::String
  local df::DataFrame

  suffix = limit == nothing ? "" : " limit $limit"
  df = ODBC.query(wrds, "$qstring$suffix")

  return df
end

wrdsQuery(wrds::ODBC.DSN, quickQuery::Symbol;
    limit::T where T<:Union{Nothing, Int} = nothing) =
    dsnQuery(wrds, QUICK_QUERIES[quickQuery], limit=limit)


#meant to be a safe, low performance way to make queries
function wrdsQuery(qstring::String;
    limit::T where T<:Union{Nothing, Int} = nothing,
    errorOnFail::Bool = true)

  local dsn::ODBC.DSN
  local df::DataFrame

  dsn = wrdsGet()

  try #try to safely make the query
    df = wrdsQuery(dsn, qstring, limit=limit)
  catch err
    wrdsDrop!(dsn) #make sure to clean up given a failure
    if !errorOnFail
      println("Query ERROR: $err")
      df = DataFrame()
    else
      error("Query ERROR: $err")
    end
  end

  wrdsDrop!(dsn)
  return df
end

wrdsQuery(quickQuery::Symbol;
    limit::T where T<:Union{Nothing, Int} = nothing) =
    wrdsQuery(QUICK_QUERIES[quickQuery], limit=limit)


#quick query examples: :libraries,  :crspTables, :crspDailyMeta, :exCrsp
function testWRDS()
  df::DataFrame = wrdsQuery(:crspaTables)

  println(df[1:200,:])

end

#@time testWRDS()
