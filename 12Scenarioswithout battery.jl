### A Pluto.jl notebook ###
# v0.17.3
################################################################################
using Markdown
using InteractiveUtils

# ╔═╡ 89ba1462-5db5-11ec-2b5a-9b7380871086
begin
	using Pkg
    PkgList = ["CSV", "DataFrames", "Gurobi", "JuMP", "Plots", "LinearAlgebra","PlutoUI","Markdown","InteractiveUtils"]
    for p in PkgList
        if Base.find_package(p) == nothing
            Pkg.add(p)
        end
	end
	using JuMP,CSV, DataFrames, Plots, PlutoUI,Markdown,InteractiveUtils, LinearAlgebra, Gurobi
end
#Paremeters#####################################################################


# ╔═╡ d3927189-5b08-483e-b678-b1c1a1c1ff95
begin
	include("readData.jl")
	BranchData, NodesData,Nodes,Pdd, Qdd, c_op = readData("E:/Msc-Skoltech/Thesis/Prof.Pozo/Meeting 2/Meeting2/Branch_data_file","E:/Msc-Skoltech/Thesis/Prof.Pozo/Meeting 2/Meeting2/Node_data_file")
end

#First Scenario (spring)

# PVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

Sbase= 1e5
#
PVDATA = CSV.read("E:/Msc-Skoltech/Thesis/Code/First-Try-withHW6/HW6_data_Demand_PV.csv",DataFrame)

PV1= 10 .*(PVDATA.PV1)/Sbase
PV2= 10 .*(PVDATA.PV2)/Sbase
PV3= 10 .*(PVDATA.PV3)/Sbase
PVplant = zeros(24,33,3)
PVplant[:,4,1] = PV1
PVplant[:,4,2] = PV2
PVplant[:,4,3] = PV3
PVplant
PV= permutedims(PVplant, (2, 1, 3))

# PVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

# PV= permutedims(PVplant, (2, 1, 3))

#
# #Second Scenario (summer)
# PV4= 14.0*(PVDATA.PV1)/100
# PV5= 14.0*(PVDATA.PV2)/100
# PV6= 14.0*(PVDATA.PV3)/100
# PVplant[:,1,4] = PV4
# PVplant[:,1,5] = PV5
# PVplant[:,1,6] = PV6
#
# #Third Scenario (fall)
# PV7= 10.0*(PVDATA.PV1)/100
# PV8= 10.0*(PVDATA.PV2)/100
# PV9= 10.0*(PVDATA.PV3)/100
# PVplant[:,1,7] = PV7
# PVplant[:,1,8] = PV8
# PVplant[:,1,9] = PV9
#
#
# #Forth Scenario (winter)
# PV10= 9.0*(PVDATA.PV1)/100
# PV11= 9.0*(PVDATA.PV2)/100
# PV12= 9.0*(PVDATA.PV3)/100
# PVplant[:,1,10] = PV10
# PVplant[:,1,11] = PV11
# PVplant[:,1,12] = PV12




Pdd= NodesData.Pd/Sbase
Qdd= NodesData.Qd/Sbase

Pmin= NodesData.Pmin/Sbase
Pmax= NodesData.Pmax/Sbase

Qmin= NodesData.Qmin/Sbase
Qmax= NodesData.Qmax/Sbase


R= BranchData.R
X= BranchData.X
NI= 33
NL= 32
# a = [0.84, 0.8, 0.772, 0.78, 0.788, 0.796, 0.788, 0.72, 0.648, 0.624, 0.56, 0.536, 0.52, 0.528, 0.536, 0.56, 0.632, 0.72, 0.88, 1.0, 1.048, 1.024, 0.96, 0.88]
#
a= [0.783041723,0.73781965,0.723014805,0.707940781,0.698788694,0.677792732,0.662718708,0.707940781,0.801076716,0.801076716,0.84333782,0.873485868,0.888559892,0.915746972,0.888559892,0.888559892,0.879407806,0.84333782,0.864333782,0.97577389,1,0.97577389,0.963930013,0.888559892]
##########for stress simulation###########################333
# a= [0.783041723,0.73781965,0.723014805,0.707940781,0.698788694,0.677792732,0.662718708,0.707940781,0.801076716,0.801076716,0.84333782,0.873485868,0.888559892,0.915746972,0.888559892,0.888559892,0.879407806,0.84333782, 1.2 .*0.864333782, 1.2 .*0.97577389, 1.0 ,0.97577389,0.963930013,0.888559892]



Pdd_ = Pdd *a'
Qdd_ = Qdd *a'
SLmax= BranchData.Smax/ 100



# *******************************************************************************
#Model

m = Model(optimizer_with_attributes(
        Gurobi.Optimizer,
        "MIPGap" => 1.7,
        "NonConvex" => 2,
        #"OutputFlag" => OutputFlag_Gurobi,
        #"max_iter" => 2000
    ))



#   ############################################################################

#Network Variables

    #define variables
    @variable(m, p[i = 1:NI,1:24,s= 1:3]); # net active withdraw at node i
    @variable(m, q[i = 1:NI,1:24,s= 1:3]); # net reactive withdraw at node i
    @variable(m, pl[l = 1:NL,1:24,s= 1:3]); #active branch power flow from i to j
    @variable(m, ql[l = 1:NL,1:24,s= 1:3]); #reactive branch power flow from i to j
    @variable(m, w[i = 1:NI,1:24,s= 1:3]); # square of voltage at node i
    @variable(m, IL[l = 1:NL,1:24,s= 1:3]); # square of current from i to j
    @variable(m, pG[i = 1:NI,1:24,s= 1:3]); #  active geneartion at node i
    @variable(m, qG[i = 1:NI, 1:24,s= 1:3]); # reactive generation at node i

	# **************************************************************************
	# @variable(m, ls[i=1:NI] >= 0) # load sheding at each node.
	# @variable(m, es[i=1:NI] >= 0) # energy spillage at each node.
	#***************************************************************************

#*******************************************************************************

# Network
@constraint( m, DistFlowEqP[l = 1:NL, t= 1:24, s= 1:3],
	pl[l,t,s] ==  p[Nodes[l, 2],t,s] + R[l] * IL[l,t,s] + sum(pl[k,t,s] for k = 1:NL if Nodes[k, 1] == Nodes[l, 2]))
@constraint(m, DistFlowEqQ[l = 1:NL,t= 1:24, s= 1:3],
	ql[l,t,s] == q[Nodes[l, 2],t,s] + X[l] * IL[l,t,s] + sum(ql[k,t,s] for k = 1:NL if Nodes[k, 1] == Nodes[l, 2]))
@constraint(m, VoltageLosses[l = 1:NL,t= 1:24, s= 1:3],
	w[Nodes[l, 2],t,s] == w[Nodes[l, 1],t,s] + (R[l]^2 + X[l]^2) * IL[l,t,s] - 2 * (R[l] * pl[l,t,s] + X[l] * ql[l,t,s]))
@constraint(m,PowerFlow_as_V_times_I[l=1:NL,t=1:24, s=1:3],(pl[l,t,s])^2+  (ql[l,t,s])^2 + ((IL[l,t,s]-w[Nodes[l,1],t,s])^2)/4 == ((IL[l,t,s]+w[Nodes[l,1],t,s])^2)/4);

# nodal net withdraw
@constraint(m, NodalEqP1[i = 1:NI, t= 1:24, s=1:3], p[i,t,s] == -pG[i,t,s] + Pdd_[i,t] - PV[i,t,s] )
@constraint(m, NodalEqQ1[i = 1:NI, t= 1:24,s=1:3], q[i,t,s] == -qG[i,t,s] + Qdd_[i,t])

# limits on load sheding and energy spillager
# @constraint(m, [i=1:NI], ls[i] <= Pdd[i])
# @constraint(m, [i=1:NI], es[i] <= pG[i])


#Technical generation limits
@constraint(m, GenLimitsP[i = 1:NI,t = 1:24,s=1:3], Pmin[i] <= pG[i,t,s] <= Pmax[i])
@constraint(m, GenLimitsQ[i = 1:NI,t= 1:24,s=1:3], Qmin[i] <= qG[i,t,s] <= Qmax[i])

#Voltage limits
@constraint(m, VoltageLimits[i = 1:NI, t=1:24, s=1:3], (0.95)^2 <= w[i,t,s] <= (1.05)^2)
#Line capacity limits
@constraint(m, FlowLimits[l = 1:NL, t= 1:24,s=1:3], (pl[l,t,s])^2 + (ql[l,t,s])^2 <= ((0.032) .^2))

#Special case for 1st node
@constraint(m, LimitsP, pG[1,1:24,1:3].==pl[1,1:24,1:3])
@constraint(m, LimitsQ, qG[1,1:24,1:3].==ql[1,1:24,1:3])
@constraint(m, SELimitsW, w[1,1:24,1:3] .== 1.00)


#Cost Vector
c= NodesData.c


COST= 1.5 .*c.*transpose(a)

#Probability of Scenarios
Prob = [0.3, 0.4, 0.3]
Prob= transpose(Prob)



# + Prob[1:3] .* 12 .* pG[1,1:24,1:3] + Prob[1:3] .* 240


# ******************************************************************************
# objective Function
# @expression(m, disp_cos1, sum(Prob[w] .*1000 .*(c[i].* (pG[i,t,w]).^2)  for i = 1:NI, t= 1:24, w= 1:3))

@expression(m, costf1, sum(c[1].* (Prob[w] .*pG[1,t,w]).^2 + Prob[w] .* 12 .* pG[1,t,w] for t= 1:24, w= 1:3) +240)
@expression(m, costf2, sum(c[18].* (Prob[w] .*pG[18,t,w]).^2 + Prob[w] .* 10.26 .* pG[18,t,w] for t= 1:24, w= 1:3) +210)
@expression(m, costf3, sum(c[22].* (Prob[w] .*pG[22,t,w]).^2 + Prob[w] .* 10.26.* pG[22,t,w] for t= 1:24, w= 1:3) +210)
@expression(m, costf4, sum(c[25].* (Prob[w] .*pG[25,t,w]).^2 + Prob[w] .* 10.26 .* pG[25,t,w] for t= 1:24, w= 1:3) +210)
@expression(m, costf5, sum(c[33].* (Prob[w] .*pG[33,t,w]).^2 + Prob[w] .* 10.26 .* pG[33,t,w] for t= 1:24, w= 1:3) +210)
@expression(m, cost_FACT, sum(1000 .* Prob[w] .*qG[18,t,w] + 1000 .* Prob[w] .*qG[33,t,w] for t= 1:24, w= 1:3) +3000)


# @expression(m, emer_cost, sum(c_ls*ls[i] + c_es*es[i] for i = 1:NI))
@objective(m, Min, costf1+costf2+costf3+costf4+costf5+ cost_FACT)

println("No model")
optimize!(m)


# ******************************************************************************
# value(disp_cost)
objective_value(m)
print(value.(w[5,1:24,1]))
pG_optimal = value.(pG)
p_optimal = value.(ql[4,18,1])
voltage_profile= sqrt.(value.(w[1:33,20,1]))
pl_optimal = value.(pl)
ql_optimal = value.(ql)

AP= (pl_optimal[1:32,20,1]).^2 + (ql_optimal[1:32,20,1]).^2

print(AP)

QG_optimal1 = 100 .*value.(qG[1,1:24,1])
QG_optimal2 = 100 .*value.(qG[18,1:24,1])
QG_optimal3 = 100 .*value.(qG[33,1:24,1])


#


p1= plot(1:33,voltage_profile,  lab = "Voltage Profile_Without ESS_ Normal",linecolor =:red4, w = 2, legend=:topleft ,xaxis= "Bus Number", yaxis= "Voltage[p.u.]", framestyle = :box, xticks =  0:1.0:33)
# p1= plot!(1:33,v_ESS_Normal,  lab = "Voltage Profile_Without ESS_Stressed",linecolor =:red1, w = 2, legend=:topleft ,xaxis= "Bus Number", yaxis= "Voltage[p.u.]", framestyle = :box, xticks =  0:1.0:33)
p1=plot!(1.05*ones(33), lab = "Max Voltage Limit = 1.05 ",line = :dash,linecolor =:darkslategrey)
p1=plot!(0.95*ones(33), lab = "Min Voltage Limit = 0.95 ",line = :dash,linecolor =:darkslategrey)










# p1=plot(1:24,Pdd_[1,1:24,1], lab = "Load", w = 3, palette = cgrad(:thermal), fill = 0, α = 0.6)
# p1= plot!(1:24,pG_optimal[1,1:24,1], lab = "Generation/ Bus1", w = 0, fillcolor=:green, fill = 0, α = 0.6, xaxis= "Time[hour]", yaxis= "P[p.u]")
p2= plot(1:24,QG_optimal1, lab = "Reactive Power_ Thermal Generator", w = 1, linecolor=:green, fillcolor=:green, fill = 0, α = 0.6, xaxis= "Time[hour]", yaxis= "Q[MVAR]",framestyle = :box, xticks =  0:1.0:24, legend=:topleft)
p2= plot!(1:24,QG_optimal2, lab = "Reactive Power_ RPC1", w = 1, fillcolor=:red, linecolor=:red, fill = 0, α = 0.5, xaxis= "Time[hour]", yaxis= "Q[MVAR]",framestyle = :box, xticks =  0:1.0:24, legend=:topleft)
p2= plot!(1:24,QG_optimal3, lab = "Reactive Power_ RPC2", w = 1, fillcolor=:blue,linecolor=:blue, fill = 0, α = 0.5, xaxis= "Time[hour]", yaxis= "Q[MVAR]",framestyle = :box, xticks =  0:1.0:24, legend=:topleft)


p1= plot(1:32, AP, lab = "Apparent Power at 20:0- without ESS ", w = 0, fillcolor=:green, fill = 0, α = 0.6, xaxis= "Line Number", yaxis= "S[MVA]", framestyle = :box, xticks =  0:1.0:33)
p1=plot!(0.001*ones(32), lab = "Max Line Capacity",line = :dash,linecolor =:red)





p1= plot(1:33,voltage_profile,  lab = "Voltage Profile_Scenario1_ESS",linecolor =:red, w = 2, legend=:topleft ,xaxis= "Time[hour]", yaxis= "Voltage[p.u.]", framestyle = :box, xticks =  0:1.0:33)

print((pl_optimal[1:32,20,1]).^2 + (ql_optimal[1:32,20,1]).^2)

print(voltage_profile)

# v_ESS_Normal = [1.0, 0.9970341993514777, 0.982972935793219, 0.9764817400004036, 0.9703610708624494, 0.9593267596277505, 0.9589326674306509,
#  0.9545569694819922, 0.952427475195657, 0.9510889688521093, 0.950639877641516, 0.95, 0.9522886442269728, 0.955325588092406, 0.9583444610303046,
#  0.961900143093166, 0.9745044857516415, 0.979439203059746,0.9965387671427081, 0.9939771446008058, 0.9937005617821777, 0.994131608709237,
#  0.9782755264361968, 0.9696525079414201, 0.9644256341731098,0.9580819429048972, 0.9565529009872775, 0.9530769256588372, 0.9511083776393392,
#  0.95, 0.9580259692047998, 0.9619413312996595, 0.9690405080292762]






# p6= plot(1:24,sqrt.(w_optimal[1,1:24,1]), lab = "V1",linecolor =:black, w = 1,ylims=(0.99,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[2,1:24,1]), lab = "V2",linecolor =:blue, w = 1,ylims=(0.99,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[3,1:24,1]), lab = "V3",linecolor =:pink, w = 1,ylims=(0.99,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[4,1:24,1]), lab = "V4",linecolor =:orange, w = 1,ylims=(0.99,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[5,1:24,1]), lab = "V5",linecolor =:purple, w = 1,ylims=(0.99,1.01),legend=:topleft)
# #


# p6= plot(1:24,sqrt.(w_optimal[1,1:24,1]), lab = "V1",linecolor =:black, w = 1,ylims=(0.995,1.01),legend=:topleft)
#
# p6= plot!(1:24,sqrt.(w_optimal[2,1:24,1]), lab = "V2",linecolor =:blue, w = 1,ylims=(0.995,1.01),legend=:topleft)
#
# p6= plot!(1:24,sqrt.(w_optimal[3,1:24,1]), lab = "V3",linecolor =:pink, w = 1,ylims=(0.995,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[4,1:24,1]), lab = "V4",linecolor =:orange, w = 1,ylims=(0.995,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[5,1:24,1]), lab = "V5",linecolor =:purple, w = 1,ylims=(0.99,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[6,1:24,1]), lab = "V6",linecolor =:brown, w = 1,ylims=(0.99,1.01),legend=:topleft)
#
# p6= plot(1:24,sqrt.(w_optimal[7,1:24,1]), lab = "V7",linecolor =:red, w = 1,ylims=(0.995,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[8,1:24,1]), lab = "V8",linecolor =:green, w = 1,ylims=(0.995,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[9,1:24,1]), lab = "V9",linecolor =:yellow, w = 1,ylims=(0.995,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[10,1:24,1]), lab = "V10",linecolor =:orange, w = 1,ylims=(0.995,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[11,1:24,1]), lab = "V11",linecolor =:black, w = 1,ylims=(0.99,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[12,1:24,1]), lab = "V12",linecolor =:purple, w = 1,ylims=(0.99,1.01),legend=:topleft)
#
# p6= plot!(1:24,sqrt.(w_optimal[13,1:24,1]), lab = "V13",linecolor =:orange, w = 1,ylims=(0.995,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[14,1:24,1]), lab = "V14",linecolor =:gray, w = 1,ylims=(0.99,1.01),legend=:topleft)
# p6= plot!(1:24,sqrt.(w_optimal[15,1:24,1]), lab = "V15",linecolor =:purple, w = 1,ylims=(0.99,1.01),legend=:topleft)
#
# p6= plot!(1:24,sqrt.(w_optimal[15,1:24,1]), lab = "V15",linecolor =:green, w = 1,ylims=(0.99,1.01),legend=:topleft)
#

savefig("E:/Msc-Skoltech/pre defence/Images/Results with ESS//apparentpower_linecaptest.png")
