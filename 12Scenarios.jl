# 12 scenarios with the Energy Storage Systems.
# Note: this is not the final version of the code implemented in the paper. 
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
	BranchData, NodesData,Nodes,Pdd, Qdd, c_op = readData("E:/.../Branch_data_file","E:/.../Node_data_file")
end

#First Scenario (spring)
Sbase= 1e5
PVDATA = CSV.read("E:/.../HW6_data_Demand_PV.csv",DataFrame)
PV1= 10 .*(PVDATA.PV1)/Sbase
PV2= 10 .*(PVDATA.PV2)/Sbase
PV3= 10 .*(PVDATA.PV3)/Sbase
PVplant = zeros(24,33,3)
PVplant[:,4,1] = PV1
PVplant[:,4,2] = PV2
PVplant[:,4,3] = PV3
PVplant
#
PVplant
PV= permutedims(PVplant, (2, 1, 3))

#Improting the data

Pdd= NodesData.Pd/Sbase
Qdd= NodesData.Qd/Sbase

Pmin= NodesData.Pmin /Sbase
Pmax= NodesData.Pmax /Sbase


Qmin= NodesData.Qmin /Sbase
Qmax= NodesData.Qmax /Sbase

SLmax= BranchData.Smax
R= BranchData.R
X= BranchData.X
NI= 33
NL= 32
# a = [0.84, 0.8, 0.772, 0.78, 0.788, 0.796, 0.788, 0.72, 0.648, 0.624, 0.56, 0.536, 0.52, 0.528, 0.536, 0.56, 0.632, 0.72, 0.88, 1.0, 1.048, 1.024, 0.96, 0.88]
a= [0.783041723,0.73781965,0.723014805,0.707940781,0.698788694,0.677792732,0.662718708,0.707940781,0.801076716,0.801076716,0.84333782,0.873485868,0.888559892,0.915746972,0.888559892,0.888559892,0.879407806,0.84333782,0.864333782,0.97577389,1,0.97577389,0.963930013,0.888559892]
##########################for stress simulation######################################
# a= [0.783041723,0.73781965,0.723014805,0.707940781,0.698788694,0.677792732,0.662718708,0.707940781,0.801076716,0.801076716,0.84333782,0.873485868,0.888559892,0.915746972,0.888559892,0.888559892,0.879407806,0.84333782, 1.2 .*0.864333782, 1.2 .*0.97577389, 1.0 ,0.97577389,0.963930013,0.888559892]

Pdd_ = Pdd *a'
Qdd_ = Qdd *a'




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


	#define variables
    @variable(m, p[i = 1:NI,1:24,s= 1:3]); # net active withdraw at node i
    @variable(m, q[i = 1:NI,1:24,s= 1:3]); # net reactive withdraw at node i
    @variable(m, pl[l = 1:NL,1:24,s= 1:3]); #active branch power flow from i to j
    @variable(m, ql[l = 1:NL,1:24,s= 1:3]); #reactive branch power flow from i to j
    @variable(m, w[i = 1:NI,1:24,s= 1:3]); # square of voltage at node i
    @variable(m, IL[l = 1:NL,1:24,s= 1:3]); # square of current from i to j
    @variable(m, pG[i = 1:NI,1:24,s= 1:3]); #  active geneartion at node i
    @variable(m, qG[i = 1:NI, 1:24,s= 1:3]); # reactive generation at node i
	@variable(m, y[1:33], Bin);
	@variable(m, xch[i = 1:NI, 1:24],Bin);
	@variable(m, xdis[i = 1:NI, 1:24],Bin);

	# **************************************************************************
	# @variable(m, ls[i=1:NI] >= 0) # load sheding at each node.
	# @variable(m, es[i=1:NI] >= 0) # energy spillage at each node.
	#***************************************************************************


#*******************************************************************************
# Battery Variables

@variable(m, Pb_ch[i=1:NI,1:24,s= 1:3])
@variable(m, Pb_dis[i=1:NI,1:24,s= 1:3])
@variable(m,Emax[i=1:NI] >= 0)
@variable(m, Qb[i=1:NI,1:24,s= 1:3])
@variable(m, e[i = 1:NI,1:24,s= 1:3])
################################################################################
# Battery Constraints

#
# @constraint(m, Pb_dismin[i=1:NI,t= 1:24,s= 1:3], 0.0 <= Pb_dis[i,t,s])
# @constraint(m, Pb_dismax[i=1:NI,t= 1:24,s= 1:3], Pb_dis[i,t,s] <= 0.8.*Emax[i])
# @constraint(m, Pbdis_binary[i=1:NI,t= 1:24,s= 1:3], Pb_dis[i,t,s] <= 0.0004 .*xdis[i,t])
# @constraint(m, Pbch_binary[i=1:NI,t= 1:24,s= 1:3], Pb_ch[i,t,s] <= 0.0004 .* xch[i,t])
# @constraint(m, Co1[i=1:NI,t= 1:24,s= 1:3],  Pb_ch[i,t,s]+Pb_dis[i,t,s] <= Emax[i])
# @constraint(m, Co123[i=1:NI,t= 1:24,s= 1:3],  Pb_ch[i,t,s]*Pb_dis[i,t,s] <= 1)
# @constraint(m,Co22[i=1:NI,t= 1:24],  xch[i,t]+xdis[i,t] <= 1 .* y[i])
# @constraint(m, batlimitQ1[i=1:NI,t= 1:24,s= 1:3], -0.0004.* y[i] <= Qb[i,t,s])
# @constraint(m, batlimitQ2[i=1:NI,t= 1:24,s= 1:3],  Qb[i,t,s] <= 0.0004 .*y[i])
# @constraint(m, maximumE_inv[i=1:NI,s= 1:3], Emax[i] <= 0.0004 .* y[i])



#
@constraint(m, Pb_dismin[i=1:NI,t= 1:24,s= 1:3], 0.0 <= Pb_dis[i,t,s])
@constraint(m, Pb_dismax[i=1:NI,t= 1:24,s= 1:3], Pb_dis[i,t,s] <= 0.8.*Emax[i])

@constraint(m, Pb_chmin[i=1:NI,t= 1:24,s= 1:3], 0.0 <= Pb_ch[i,t,s])
@constraint(m, Pb_chmax[i=1:NI,t= 1:24,s= 1:3], Pb_ch[i,t,s] <=  Emax[i])

@constraint(m, Pbdis_binary[i=1:NI,t= 1:24,s= 1:3], Pb_dis[i,t,s] <= 2000 .*xdis[i,t])
@constraint(m, Pbch_binary[i=1:NI,t= 1:24,s= 1:3], Pb_ch[i,t,s] <= 2000 .* xch[i,t])
@constraint(m, Co1[i=1:NI,t= 1:24,s= 1:3],  Pb_ch[i,t,s]+Pb_dis[i,t,s] <= Emax[i])
@constraint(m, Co123[i=1:NI,t= 1:24,s= 1:3],  Pb_ch[i,t,s] .*Pb_dis[i,t,s] <= 1)
@constraint(m,Co22[i=1:NI,t= 1:24],  xch[i,t] +xdis[i,t] <= 1200 .* y[i])
@constraint(m,Co222[i=1:NI,t= 1:24],  xch[i,t] + xdis[i,t] <= 1)
@constraint(m, batlimitQ1[i=1:NI,t= 1:24,s= 1:3], -1200 .* y[i] <= Qb[i,t,s])
@constraint(m, batlimitQ2[i=1:NI,t= 1:24,s= 1:3],  Qb[i,t,s] <= 1200 .*y[i])
@constraint(m, maximumE_inv[i=1:NI,s= 1:3], Emax[i] <= 1200 .* y[i])





# @constraint(m, Pb_dismin[i=1:NI,t= 1:24, s= 1:3], 0.0 <= Pb_dis[i,t,s])
# @constraint(m, Pb_dismax[i=1:NI,t= 1:24, s= 1:3], Pb_dis[i,t,s] <= 0.8.*Emax[i])
#
# @constraint(m, Pb_chmin[i=1:NI,t= 1:24, s=1:3], 0.0 <= Pb_ch[i,t,s])
# @constraint(m, Pb_chmax[i=1:NI,t= 1:24, s=1:3], Pb_ch[i,t,s] <=  Emax[i])
#
# @constraint(m, Pbdis_binary[i=1:NI,t= 1:24, s=1:3], Pb_dis[i,t,s] <= 2000 .*xdis[i,t])
# @constraint(m, Pbch_binary[i=1:NI,t= 1:24, s=1:3], Pb_ch[i,t,s] <= 2000 .* xch[i,t])
# @constraint(m, Co1[i=1:NI,t= 1:24, s=1:3],  Pb_ch[i,t,s]+Pb_dis[i,t,s] <= Emax[i])
# @constraint(m, Co123[i=1:NI,t= 1:24, s=1:3],  Pb_ch[i,t,s]*Pb_dis[i,t,s] <= 1)
# @constraint(m,Co22[i=1:NI,t= 1:24],  xch[i,t]+xdis[i,t] <= 2000 .* y[i])
# @constraint(m, batlimitQ1[i=1:NI,t= 1:24, s=1:3], -1200 .* y[i] <= Qb[i,t,s])
# @constraint(m, batlimitQ2[i=1:NI,t= 1:24, s=1:3],  Qb[i,t,s] <=1200 .*y[i])
# @constraint(m, maximumE_inv[i=1:NI, s=1:3], Emax[i] <= 12000 .* y[i])



################################################################################
@constraint(m,sum(y[1:33]) <= 1)  # no more than 1 battery
@constraint(m,Socmin[i=1:NI,t= 1:24, s=1:3], 0.0 <= e[i,t,s])
@constraint(m,Socmax[i = 1:NI,t=1:24, s=1:3], e[i,t,s] <= 100)

@constraint(m,SOC1[i=1:NI,t= 1:23, s=1:3], e[i,t.+1,s] == e[i,t,s] - 1.11 *Pb_dis[i,t.+1,s] + 0.85 *Pb_ch[i,t.+1,s])
@constraint(m, efirststatus[i=1:NI,t= 1:24, s=1:3], e[i,1,s] == e[i,24,s] - 1.11 *Pb_dis[i,24,s] + 0.85 *Pb_ch[i,24,s])
# @constraint(m,ce[i = 1:NI,t=1:24, s= 1:3], e[i,t,s] == 0.9 .* Emax[i]);
#*******************************************************************************
# Network Constraints
@constraint( m, DistFlowEqP[l = 1:NL, t= 1:24, s=1:3],
	pl[l,t,s] ==  p[Nodes[l, 2],t,s] + R[l] * IL[l,t,s] + sum(pl[k,t,s] for k = 1:NL if Nodes[k, 1] == Nodes[l, 2]))
@constraint(m, DistFlowEqQ[l = 1:NL,t= 1:24, s=1:3],
	ql[l,t,s] == q[Nodes[l, 2],t,s] + X[l] * IL[l,t,s] + sum(ql[k,t,s] for k = 1:NL if Nodes[k, 1] == Nodes[l, 2]))
@constraint(m, VoltageLosses[l = 1:NL,t= 1:24, s=1:3],
	w[Nodes[l, 2],t,s] == w[Nodes[l, 1],t,s] + (R[l]^2 + X[l]^2) * IL[l,t,s] - 2 * (R[l] * pl[l,t,s] + X[l] * ql[l,t,s]))
@constraint(m,PowerFlow_as_V_times_I[l=1:NL,t=1:24,s= 1:3],(pl[l,t,s])^2+  (ql[l,t,s])^2 + ((IL[l,t,s]-w[Nodes[l,1],t,s])^2)/4 == ((IL[l,t,s]+w[Nodes[l,1],t,s])^2)/4);




# nodal net withdraw (demand which is varying by Solar generation)
@constraint(m, NodalEqP1[i = 1:NI, t= 1:24, s=1:3], p[i,t,s] == -pG[i,t,s] + Pdd_[i,t] + Pb_ch[i,t,s]- Pb_dis[i,t,s] - PV[i,t,s])
@constraint(m, NodalEqQ1[i = 1:NI, t= 1:24, s=1:3], q[i,t,s] == -qG[i,t,s] + Qdd_[i,t]- Qb[i,t,s])

# limits on load sheding and energy spillager
# @constraint(m, [i=1:NI], ls[i] <= Pdd[i])
# @constraint(m, [i=1:NI], es[i] <= pG[i])


#Technical generation limits
@constraint(m, GenLimitsP[i = 1:NI,t = 1:24, s=1:3], Pmin[i] <= pG[i,t,s] <= Pmax[i])
@constraint(m, GenLimitsQ[i = 1:NI,t= 1:24, s=1:3], Qmin[i] <= qG[i,t,s] <= Qmax[i])

#Voltage limits
@constraint(m, VoltageLimits[i = 1:NI, t=1:24, s=1:3], (0.95)^2 <= w[i,t,s] <= (1.05)^2)

Line_factors = [1,1,1,0.2]
Smax = Line_factors .* SLmax[1:4]
# @constraint(m, FlowLimits2[l = 1:NL, t= 1:24,s=1:3], (pl[l,t,s])^2 + (ql[l,t,s])^2 <= (Smax[l]).^2)
@constraint(m, FlowLimits2[l = 1:NL, t= 1:24, s=1:3], (pl[l,t,s])^2 + (ql[l,t,s])^2 <= (1).^2)



#Special case for 1st node
@constraint(m, LimitsP, pG[1,1:24,1:3].==pl[1,1:24,1:3])
@constraint(m, LimitsQ, qG[1,1:24,1:3].==ql[1,1:24,1:3])
@constraint(m, SELimitsW, w[1,1:24,1:3] .== 1.00)


#COST

Prob = [0.3, 0.4, 0.3]
Prob= transpose(Prob)
# b_cost= [0.5,0.5,0.5,0.5,0.5]
# ai_= [561,0.0, 310,310,78]
# ai= transpose(a).*ai_
# aii= ai[:,:,:]
#
# bi= [7.29,0.0, 7.85,7.85,7.97]
# bi= transpose(a).*bi
# bii=bi[:,:,:]
# ci_=[0.00156,0.0, 0.00194,0.00194,0.00482]
# ci= transpose(a) .*c
# cii=ci[:,:,:]
# Prob = [0.12, 0.049, 0.089,0.074,0.083,0.066,0.046,0.097,0.112,0.08,0.089,0.095]


# c
# c= c[:,:]

#
# ci[1:15,1:24] .*pG[1:15,1:24,1:3]
# ******************************************************************************

#Cost

c= NodesData.c
c =c
cost_coefficients = [0.035,0.028,0.029,0.042,0.03,0.04,0.07,0.09,0.09,0.06,0.06,0.05,0.04,0.029,0.026,0.028,0.0319,0.024,0.025,0.024,0.024,0.024,0.024,0.024]
COST= c.*transpose(cost_coefficients)

#Objective function

@expression(m, costf1, sum(0.3.* (Prob[w] .*pG[1,t,w]).^2 + Prob[w] .* 1200 .* pG[1,t,w] for t= 1:24, w= 1:3) +240)
# @expression(m, costf2, sum(c[18].* (Prob[w] .*pG[18,t,w]).^2 + Prob[w] .* 10.26 .* pG[18,t,w] for t= 1:24, w= 1:3) +210)
# @expression(m, costf3, sum(c[22].* (Prob[w] .*pG[22,t,w]).^2 + Prob[w] .* 10.26.* pG[22,t,w] for t= 1:24, w= 1:3) +210)
# @expression(m, costf4, sum(c[25].* (Prob[w] .*pG[25,t,w]).^2 + Prob[w] .* 10.26 .* pG[25,t,w] for t= 1:24, w= 1:3) +210)
# @expression(m, costf5, sum(c[33].* (Prob[w] .*pG[33,t,w]).^2 + Prob[w] .* 10.26 .* pG[33,t,w] for t= 1:24, w= 1:3) +210)

@expression(m, operation_ESS, sum(0.5 .* Pb_ch[i,t,w] + 0.5 .*Pb_dis[i,t,w]  for i = 1:NI, t= 1:24, w= 1:3))
@expression(m, investment_ESS, sum(25000 .* Emax[i] for i = 1:NI))
# @expression(m, costf_ESS, (sum(0.001 .* (Pb_ch[i,t,w]) + 0.0025 .* Emax[i]  for i = 1:NI, t= 1:24, w= 1:3)))

# @expression(m, costf_ESS,  sum(0.005 .*(0.1 .* (Pb_ch[i,t,w]).^2 + 0.25 .* (Emax[i]).^2  for i = 1:NI, t= 1:24, w= 1:3)))
# @expression(m, costf_ESS,  sum(0.005 .*(0.1 .* (Prob[w] .* (Pb_ch[i,t,w]) + 0.25 .* Emax[i]  for i = 1:NI, t= 1:24, w= 1:3)))
# @expression(m, costf_ESS2, sum((250 .* (Emax[i]))  for i = 1:NI))
# @expression(m, emission, 23 .* (sum(3.21 .* (0.00004.* (Prob[w] .*pG[1,t,w]).^2 + Prob[w] .* 0.3 .* pG[1,t,w]) for t= 1:24, w= 1:3) +45))
# @expression(m, cost_Q, sum(1000 .* Prob[w] .*Qb[i,t,w] for  i = 1:NI, t= 1:24, w= 1:3))


@objective(m, Min, costf1+operation_ESS+investment_ESS)


# @objective(m, Min, sum(Prob[w] .*1000 .*(c[i].* (pG[i,t,w]).^2) + Prob[w] .* 1000 .*(Pb_ch[i,t,w]).^2  + 0.0 .* (Emax[i]) for i in 1:NI, t in 1:24, w in 1:3))




println("No model")
optimize!(m)
# ******************************************************************************
# value(disp_cost)
ov = objective_value(m)
w_optimal = value.(w)
pG_optimal = 100 .*value.(pG)
QG_optimal = 100 .*value.(qG[1,1:24,1])
Qb_optimal = 100 .*value.(Qb[28,1:24,1])
pl_optimal = value.(pl)
ql_optimal = value.(ql)
qG_optimal = value.(pl)
qG_optimal = value.(ql)

Pb_optimal1 = value.(Pb_ch)
Pb_optimal2 = value.(Pb_dis)

AP= (pl_optimal[1:32,20,1]).^2 + (ql_optimal[1:32,20,1]).^2

demand1= value.(p[4,1:24,1])
demand2= value.(p[4,1:24,2])
demand3= value.(p[4,1:24,3])


Emax_optimal = value.(Emax)


(value.(pl[1:14,21,1])).^2+(value.(ql[1:14,21,1])).^2
voltage_profile= sqrt.(value.(w[1:33,20,1]))
voltage_profile2= sqrt.(value.(w[1:33,20,2]))
voltage_profile3= sqrt.(value.(w[1:33,20,3]))
# voltage_profile2= sqrt.(value.(w[1:33,20,1]))




############################## Plots ######################################

# # ╔═╡ 84ee4bcc-7a61-4710-9440-b1627759080b

#Voltage Profile
p1= plot(1:33,voltage_profile,  lab = "Voltage Profile_With ESS_ Normal",linecolor =:blue4, w = 2, legend=:topleft ,xaxis= "Bus Number", yaxis= "Voltage[p.u.]", framestyle = :box, xticks =  0:1.0:33)
p1= plot!(1:33,transpose(V_ESS_Normala),  lab = "Voltage Profile_With ESS_Stressed",linecolor =:deepskyblue1, w = 2, legend=:topleft ,xaxis= "Bus Number", yaxis= "Voltage[p.u.]", framestyle = :box, xticks =  0:1.0:33)
p1=plot!(1.05*ones(33), lab = "Max Voltage Limit = 1.05 ",line = :dash,linecolor =:darkslategrey)
p1=plot!(0.95*ones(33), lab = "Min Voltage Limit = 0.95 ",line = :dash,linecolor =:darkslategrey)




#Generation
# p1= plot(1:24,pG_optimal[1,1:24,1], lab = "Generation/ Bus1", w = 1, fillcolor=:green, fill = 0, α = 0.6, xaxis= "Time[hour]", yaxis= "P[MW]",framestyle = :box)
# p1= plot(1:24,pG_optimal[1,1:24,1], lab = "Generation/ Bus1", w = 1, fillcolor=:green, fill = 0, α = 0.6, xaxis= "Time[hour]", yaxis= "P[MW]")
#Generation_reactive power
# q1= plot(1:24,Qb_optimal, lab = "ESS_Reactive_Power", w = 1, fillcolor=:seagreen1, fill = 0, α = 0.8, xaxis= "Time[hour]", yaxis= "P[MW]")
# q1= plot!(1:24,QG_optimal, lab = "Thermal_Gen_Reactive_Power", w = 1, fillcolor=:darkorchid4, linecolor=:darkorchid4, fill = 0, α = 0.8,legend=:topleft, xaxis= "Time[hour]", yaxis= "Q[MVar]", xticks =  0:1.0:24)

#Load Profile
# p1= plot!(1:24,3.17 .* a, lab = "Normalized Load Profile",linecolor =:teal, w = 1,legend=:topleft,xaxis= "Time[hour]", yaxis= "" , framestyle = :box,xticks =  0:1.0:24)

#ESS
p1= plot(2:23,(-100) .*Pb_optimal1[27, 2:23,1], lab = "Charging Power",w = 2, fillcolor=:grey31,linecolor=:grey31, fill = 0, α = 0.8,xaxis= "Time[hour]", yaxis= "P[MW]",xticks =  0:2.0:24)
p1= plot!(2:23,100 .* Pb_optimal2[27,2:23,1], lab = "Discharging Power",w = 2, fillcolor=:teal,linecolor=:teal, fill = 0, α = 0.8,legend=:top, framestyle = :box,xticks =  0:2.0:24)

#Apparent power with ESS
p1= plot(1:32, AP, lab = "Apparent Power at 20:00,Stressed Condition,without ESS", w = 0, fillcolor=:royalblue1, fill = 0, α = 0.6, xaxis= "Line Number", yaxis= "S[MVA]/100", framestyle = :box, xticks =  0:1.0:33)
p1= plot!(1:32, AP_withoutESS, lab = "Apparent Power at 20:00,Stressed Condition,ESS", w = 0, fillcolor=:navyblue, fill = 0, α = 0.6, xaxis= "Line Number", yaxis= "S[MVA]/100", framestyle = :box, xticks =  0:1.0:33)
p1=plot!(0.002*ones(32), lab = "Max Line Capacity",line = :dash,linecolor =:red)



#DEMAND
d3= plot(1:24,100 .*demand3, lab = "Demand/Scenario-1/Bus No.4", w = 1.2, linecolor=:green, xaxis= "Time[hour]", yaxis= "P[MW]",legend=:bottomleft, framestyle = :box,xticks =  0:1.0:24)
d2= plot!(1:24,100 .*demand2, lab = "Demand/Scenario-2/Bus No.4", w = 1.2, linecolor=:blue,xaxis= "Time[hour]", yaxis= "P[MW]",framestyle = :box,xticks =  0:1.0:24)
d1= plot!(1:24,100 .*demand1, lab = "Demand/Scenario-3/Bus No.4", w = 1.2, linecolor=:red, xaxis= "Time[hour]", yaxis= "P[MW]",framestyle = :box,xticks =  0:1.0:24)
# d2= plot!(1:24,demand2, lab = "Generation/ Bus1", w = 1, fillcolor=:blue, fill = 0, α = 0.6, xaxis= "Time[hour]", yaxis= "P[MW]",framestyle = :box)
# d3= plot!(1:24,demand3, lab = "Generation/ Bus1", w = 1, fillcolor=:green, fill = 0, α = 0.6, xaxis= "Time[hour]", yaxis= "P[MW]",framestyle = :box)





print((pl_optimal[1:32,20,1]).^2 + (ql_optimal[1:32,20,1]).^2)

savefig("E:/Msc-Skoltech/pre defence/Images/Results with ESS//apparentpowerstressedcondition.png")


# V_without_ESS= [1,0.996222662,0.977808921,0.968085262,0.958580136,0.935036648,0.930328486,0.924350078,0.917154184,0.910710635,0.909821007,0.908347407,0.902655486,0.900730804,0.900151039,0.9,0.900793973,0.901838402,0.995726914,0.993164037,0.992887361,0.993318759,0.973086023,0.964415996,0.959160442,0.932652969,0.9295337,0.915824748,0.906141839,0.902199149,0.9,0.900084225,0.90163928]






# PV1= plot(1:24,PVplant[1:24,1,3], lab = "PV1-Sunny day", w = 1, fillcolor=:yellow, linecolor=:yellow, fill = 0, α = 1)
# PV2= plot!(1:24,PVplant[1:24,1,2], lab = "PV2-Cloudy day", w = 1, fillcolor=:orange,linecolor=:orange, fill = 0, α = 1)
# PV3= plot!(1:24,PVplant[1:24,1,1], lab = "PV3-Rainy day", w = 1, fillcolor=:red, linecolor=:red,fill = 0, α = 1)


