#############################
# Reading Power flow data from matpower (.m) file
#############################


    bus = ["Bus No", "type", "Pd", "Qd", "Gs", "Bs", "area", "Vm", "Va", "baseKv", "zone", "Vmax", "Vmin"]
    gen = ["bus", "Pg", "Qg", "Qmax", "Qmin", "Vg", "mBase", "status", "Pmax", "Pmin", "Pc1", "Pc2", "Qc1min", "Qc1max",
        "Qc2min", "Qc2max", "ramp_agc", "ramp_10", "ramp_30", "ramp_q", "apf"]
    branch = ["fbus", "tbus", "r", "x", "b", "rateA", "rateB", "rateC", "ratio", "angle", "status", "angmin", "angmax"]

function ReadMatData(file::String)
  if endswith(file, ".m")
    NetworkData = parse_data(file)
  end
end

function Str_to_Var(s::String)
   s = Symbol(s)
   @eval (($s))
 end

function parse_data(datafile :: String)
  global Network = Dict()
  datalines = readlines(datafile)
  Index = 1
  while Index <= length(datalines)
    if length(strip(datalines[Index]))<=0 || strip(datalines[Index])[1] == '%'
      Index = Index +1
      continue
    end
    if contains(datalines[Index], "function mpc")
      Network["Name"] = strip(split(datalines[Index], "=")[2], ['\n', ';', '"'])
    elseif contains(datalines[Index], "mpc.baseMVA")
      Network["Base MVA"]= strip(split(datalines[Index], "=")[2], ['\n', ';', '"'])
      elseif contains(datalines[Index], "mpc.") && contains(datalines[Index], "[")
            if contains(datalines[Index], "bus") || contains(datalines[Index], "gen") || contains(datalines[Index], "branch")
                if !contains(datalines[Index], "gencost")
                  SParameterType = strip(replace(split(datalines[Index])[1], "mpc.", ""))
                  ParameterType = Str_to_Var(SParameterType)
                  DataValue, Index = parse_MatrixData(datalines, Index)
                  MakeDict(DataValue,ParameterType,SParameterType)
                end
            end
      end
    Index=Index + 1
  end
end

function parse_MatrixData(datalines::Array{String,1}, Index::Int)
    MatrixData = Array[];
    while Index <= length(datalines)
      if length(strip(datalines[Index]))<=0 || contains(datalines[Index],"%") || contains(datalines[Index],"[")
        Index = Index + 1
        continue
      end
      contains(datalines[Index],"]") ? break : ""
        push!(MatrixData, split(replace(datalines[Index], ";","")))
      Index = Index + 1
    end
    return MatrixData, Index
end

function MakeDict(ParameterValue, ParameterType, SParameterType)
    Parameters=Dict();
    for irow in 1:length(ParameterValue)
        TypeDict = Dict()
        for icol in 1:length(ParameterType)
         TypeDict[ParameterType[icol]]= parse(Float64, ParameterValue[irow][icol])
        end
        Parameters["$irow"]=TypeDict
    end
    Network[SParameterType]=Parameters
end

ReadMatData("Case30Bus.m")

Network
