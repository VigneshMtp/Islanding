#Islanding
using JuMP
using Gurobi
using PowerModels
global data=PowerModels.parse_file("Case30Bus.m")

IslandingOutFile = open("OuputFile.txt", "w");close(IslandingOutFile)

T_edge = Array(Vector{Int64}, length(data["bus"])) # Set of nodes connected to particular bus
D_edge1 = Array(Vector{Int64}, length(data["bus"]))
D_edge2 = Array(Vector{Int64}, length(data["bus"]))
T_x = Array(Vector{Float64}, length(data["bus"])) # Reactance of edges connected to node

[T_edge[j]=[] for j in 1:length(data["bus"])]
[T_x[j]=[] for j in 1:length(data["bus"])]

[D_edge1[j]=[] for j in 1:length(data["bus"])]
[D_edge2[j]=[] for j in 1:length(data["bus"])]

[push!(T_edge[data["branch"]["$i"]["f_bus"]], data["branch"]["$i"]["t_bus"]) for i in 1:length(data["branch"])]
[push!(T_edge[data["branch"]["$i"]["t_bus"]], data["branch"]["$i"]["f_bus"]) for i in 1:length(data["branch"])]

[push!(D_edge1[data["branch"]["$i"]["f_bus"]], data["branch"]["$i"]["t_bus"]) for i in 1:length(data["branch"])]
[push!(D_edge2[data["branch"]["$i"]["t_bus"]], data["branch"]["$i"]["f_bus"]) for i in 1:length(data["branch"])]

[push!(T_x[data["branch"]["$i"]["f_bus"]], data["branch"]["$i"]["br_x"]) for i in 1:length(data["branch"])]
[push!(T_x[data["branch"]["$i"]["t_bus"]], data["branch"]["$i"]["br_x"]) for i in 1:length(data["branch"])]

alp_val = 100

for N_isl in 2 : 6
  for alp_ind in 1 : 11
    IslandingOutFile = open("OuputFile.txt", "a")
    DCOPF = Model(solver = GurobiSolver(MIPGap = 1e-2, TimeLimit = 800, OutputFlag=0))
    @defVar(DCOPF, pgen[1:length(data["bus"])]>=0)
    @defVar(DCOPF, -4 <= delta[1:length(data["bus"])]<=4)
    @defVar(DCOPF, x[1:length(data["bus"]), 1:N_isl], Bin)
    @defVar(DCOPF, y[1:length(data["bus"]), 1:length(data["bus"])], Bin)
    @defVar(DCOPF, r[1:length(data["bus"])]>=0)
    @defVar(DCOPF, w[1:length(data["bus"])+1, 1:length(data["bus"]), k in 1:N_isl]>=0)
    @defVar(DCOPF, sg[1:length(data["bus"]), k in 1:N_isl], Bin)
    @defVar(DCOPF, beta[1:length(data["bus"]), 1:length(data["bus"]), k in 1:N_isl], Bin)
    @defVar(DCOPF, gamma[1:length(data["bus"]), 1:length(data["bus"]), k in 1:N_isl], Bin)
    @defVar(DCOPF, Size>=0)
    @setObjective(DCOPF, Min, sum(10*r[i]*100 for i in 1:length(data["bus"]))+(alp_ind-1)*alp_val*(Size)/10)

    #Power Flow Constraints
    @addConstraint(DCOPF, CheckConst[j in 1:length(data["bus"])], pgen[j]*(data["bus"]["$j"]["bus_type"]-2)>=0)
    @addConstraint(DCOPF, PBconst[i in 1:length(data["bus"])], pgen[i] -
        sum(1/T_x[i][j]*(delta[i] - delta[T_edge[i][j]])y[i,T_edge[i][j]] for j in 1:length(T_edge[i])) - (data["bus"]["$i"]["pd"] - r[i])==0)
    @addConstraint(DCOPF, LoadConst[i in 1:length(data["bus"])], (data["bus"]["$i"]["pd"] - r[i])>=0)
    @addConstraint(DCOPF, delconst1[j in 1:length(data["bus"])], delta[j]<=4*(1-sum(sg[j,k] for k in 1:N_isl)))
    @addConstraint(DCOPF, delconst2[j in 1:length(data["bus"])], delta[j]>=-4*(1-sum(sg[j,k] for k in 1:N_isl)))
    @addConstraint(DCOPF, GenConst1[i in 1:length(data["gen"])],pgen[data["gen"]["$i"]["gen_bus"]]>=data["gen"]["$i"]["pmin"])
    @addConstraint(DCOPF, GenConst2[i in 1:length(data["gen"])],pgen[data["gen"]["$i"]["gen_bus"]] <= data["gen"]["$i"]["pmax"])

    #Islanding Constraints
    @addConstraint(DCOPF, IslandConst1[i in 1:length(data["bus"])], sum(x[i,k] for k in 1:N_isl)==1)
    @addConstraint(DCOPF, IslandConst21[l in 1:length(data["branch"])], y[data["branch"]["$l"]["f_bus"],data["branch"]["$l"]["t_bus"]] ==
        sum(x[data["branch"]["$l"]["f_bus"],k]*x[data["branch"]["$l"]["t_bus"],k] for k in 1:N_isl))
    @addConstraint(DCOPF, IslandConst22[l in 1:length(data["branch"])], y[data["branch"]["$l"]["t_bus"],data["branch"]["$l"]["f_bus"]] ==
        sum(x[data["branch"]["$l"]["t_bus"],k]*x[data["branch"]["$l"]["f_bus"],k] for k in 1:N_isl))
    @addConstraint(DCOPF, SizeConst[k in 1:N_isl], Size>=sum(x[i,k] for i in 1:length(data["bus"]))+sum(sg[j,k] for j in 1:length(data["bus"]), k in 1:N_isl))

    # Network Connectivity
    @addConstraint(DCOPF, NCL_Const[j in 1:length(data["bus"]), k in 1:N_isl], w[length(data["bus"])+1,j,k]==sum(beta[jj,j,k] for jj in 1:length(data["bus"])))
    @addConstraint(DCOPF, NCL_Const1[jj in 1:length(data["bus"]), j in 1:length(data["bus"]), k in 1:N_isl], beta[jj,j,k]<=x[jj,k])
    @addConstraint(DCOPF, NCL_Const3[jj in 1:length(data["bus"]), j in 1:length(data["bus"]), k in 1:N_isl], beta[jj,j,k]<=sg[j,k])
    @addConstraint(DCOPF, NCL_Const4[jj in 1:length(data["bus"]), j in 1:length(data["bus"]), k in 1:N_isl], beta[jj,j,k]>=x[jj,k]+sg[j,k]-1)
    @addConstraint(DCOPF, NCL_Const5[j in 1:length(data["bus"]), k in 1:N_isl], sg[j,k]<=x[j,k])
    @addConstraint(DCOPF, NC_Const2[k in 1:N_isl], sum(sg[j,k] for j in 1:length(data["bus"]))<=1)
    @addConstraint(DCOPF, NC_Const3[j in 1:length(data["bus"]), k in 1:N_isl], x[j,k]
       == sum(w[D_edge2[j][i],j,k] for i in 1:length(D_edge2[j]))-sum(w[j,D_edge1[j][i],k] for i in 1:length(D_edge1[j]))
       + w[length(data["bus"])+1,j,k])
    @addConstraint(DCOPF, NClim_Const1[i in 1:length(data["bus"]), j in 1:length(data["bus"]), k in 1:N_isl],w[i,j,k]<=length(data["bus"])*x[i,k])
    @addConstraint(DCOPF, NClim_Const2[i in 1:length(data["bus"]), j in 1:length(data["bus"]), k in 1:N_isl],w[i,j,k]<=length(data["bus"])*x[j,k])

    solve(DCOPF, suppress_warnings =true)

    print(IslandingOutFile, N_isl, "\t",(alp_ind-1)*alp_val/10, "\t" ) # Maximum island bounds & Severity level

    # Subnetwork nodes
    for k in 1:N_isl
      for i in 1:length(data["bus"])
        round((getvalue(x[i,k])))==1 ? print(IslandingOutFile, i,","):""
      end
      print(IslandingOutFile, "\n\t\t")
    end
    print(IslandingOutFile, "\t\t")
    print(IslandingOutFile, round(sum(getvalue(r[:]))*100),"\t") # Load shedding

    # Generators' output
    for i in 1:length(data["bus"])
      data["bus"]["$i"]["bus_type"]>=2 ? print(IslandingOutFile,round(getvalue(pgen[i])*100), ","):""
    end
    print(IslandingOutFile, "\n")
    close(IslandingOutFile)
  end
end
