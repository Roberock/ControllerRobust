using Pkg
import Pkg;  
Pkg.add("StatsPlots")
Pkg.add("ControlSystems"); 
Pkg.add("Distributed");
Pkg.add("StatsBase");
Pkg.add("Distributions");
Pkg.add("CSV")
# Pkg.add("ControlToolbox")

 
using ControlSystems, StatsPlots
using Distributed, StatsBase, Distributions 
using CSV, DataFrames

include(string(pwd(), "\\src\\build_plant.jl"));
include(string(pwd(), "\\src\\g_controller.jl"));
include(string(pwd(), "\\src\\tf2ss.jl"));
include(string(pwd(), "\\src\\MonteCarloController.jl"));
include(string(pwd(), "\\src\\EstimateReliabilityScores.jl"));
include(string(pwd(), "\\data\\MyDGM.jl"));## my Data-Generating-MechanismDGM
#####################################################################################
####################  Define contorl problem and run example ######################## 
#####################################################################################

## Parameter δ (uncertain parameters)
δref=[1, 1, 1, 0, 1,0.001]      ## example: a reference value for the uncertain parameters δ
δnames=["mass1" "mass2" "kl" "kn" "λ" "τ" ]     ## names of the uncertain factors
# admissible ranges
δlims=[0.1    2;        0.1    2;
       0.05   1.75;     -1     1;
       0.2    1.8;      0.0001 0.3];             

# nominal design      
dnom=[-0.1324, 0.3533, 0.6005, 0.0728, 0.5503, 1.4175, 2.6531, 2.4802, 1];
d_names =["a1" "a2" "a3" "a4" "a5" "b1" "b2" "b3" "b4"]     ## names of design parameters 

### Example: build the system matrices (state-space representation)
A, B, C, D = build_plant(δref); 

### Example: define the controller convering from trasfer function to system state matrices  
Ac, Bc, Cc, Dc=tf2ss(dnom[1:4],dnom[5:end]);  

### Example: call solver and compute reliability function and system response 
tvec=0:0.01:25;
g1_nom,g2_nom,g3_nom, y_nominal =g_controller(dnom,δref,tvec,true);
 
## Example: Monte Carlo simulation  (for random samples within δlims)
Ns=1500
δuniform =  δlims[:,1] .+ (δlims[:,2]- δlims[:,1]).*rand(6,Ns);
# δgaussian=  MyDGM(δref, δlims, Ns) 

Gref ,Yref =  MonteCarloController(dnom,Array(δuniform'),tvec)

 ## Exanoke plot system response   
plot(tvec,  Yref' , alpha=0.5 ,legend = false, xlabel="time", ylabel="x2(t)") 
plot(tvec[1:500],  Yref[:,1:500]' , alpha=0.5 ,legend = false, xlabel="time", ylabel="x2(t)") 

## Example:  reliability scores from given-data
Pf_ind, Pf_all, Sev_ind, Sev_all, W, IsF_all, IsF_ind= EstimateReliabilityScores(Gref);
 

######################################################################
####################  Analyze one of 4 designs ####################### 
######################################################################
# Candidate controller design    
#   𝒮𝒫1(𝒟,0) 𝒮𝒫1(𝒟,.05) 𝒮𝒫2(𝒟) Nominal
# 𝑎4  0.2238   0.5375    0.7600  0.5503
# 𝑎3  0.6811   1.3346    1.9491  1.4175
# 𝑎2  3.1275   2.4206    3.0497  2.6531
# 𝑎1  2.3615   2.1689    2.7344  2.4802
# 𝑎0  1.1833   0.8084    1.0594  1.0000
# 𝑏3 −0.0982   2.4802   −0.0831 -0.1324  
# 𝑏2  0.4702   0.6146    0.6358  0.3533
# 𝑏1  0.5886   0.5265    0.7752  0.6005
# 𝑏0  0.0777   0.0716    0.0981  0.0728
 
Design_Selection=1
if Design_Selection==1 # same as the nominal
    d_new=[-0.1324, 0.3533, 0.6005, 0.0728, 0.5503, 1.4175, 2.6531, 2.4802, 1];  
elseif Design_Selection==2 # desing from min{ \hat{F}_w^{-1}(1) }
    d_new=[ -0.103539257684432,	0.452027845175231,	0.553521292381580,	0.0731453358298691,	0.285291345952459, 0.703095280659925,	3.45687805905316,	2.49234017067302,	1.16716573966267]; # 0.717520994282477
elseif Design_Selection==3 # from min { \hat{F}_w^{-1}(0.95) }
    d_new=[-0.0276450664037319, 0.613233915821064, 0.528302902131644, 0.0732150560642133, 0.585865170685746, 1.35042413219186, 2.51075569540054, 2.17663965704309, 0.827018013978000];
elseif Design_Selection==4 # from  min { Pf }
    d_new=[-0.0831455337335365,	0.635784535253079,	0.775222237017883,	0.0981477664191349,	0.760033487719863,	1.94907797910329,	3.04966829435795,	2.73439891631291,	1.05935063251108]; #  gamma: 0.2218
end
 
## load data (samples of the uncertain factors)

DataSetType =1
Nuber_of_samples=10^3
if DataSetType==1
    File_path=string(pwd(), "\\data\\Theta_small.csv"); ## small uncertainty
    df = CSV.read(File_path,DataFrame);
    δ_samples=Matrix(df);
elseif DataSetType==2
    File_path=string(pwd(), "\\data\\Theta_large.csv"); ## large uncertainty
    df = CSV.read(File_path,DataFrame);
    δ_samples=Matrix(df);
#elseif DataSetType==3 
#    δ_samples=MyDGM(Nuber_of_samples); % your own data generating mechanism
end

 
## Monte Carlo simulation for d_new, failure probability, and severity  
G, Y_samples  =  MonteCarloController(d_new,δ_samples,tvec)

Pf_ind, Pf_all, Sev_ind, Sev_all, W, IsFail_all, IsF_ind= EstimateReliabilityScores(G);  

Y_safe=Y_samples[vec(IsFail_all.==0),:]';
Y_fail=Y_samples[vec(IsFail_all),:]';
 ## plot system response   
plot(tvec,  Y_safe , alpha=0.1 ,legend = false, xlabel="time", ylabel="x2(t)")
plot!(tvec, Y_fail , color =[:red], legend = false, xlabel="time", ylabel="x2(t)")

## plot samples of g1,g2,g3   
 l = @layout [a ; b c]

 p1 = scatter( G[vec(IsFail_all),1], G[vec(IsFail_all),2],labels="fail" ,markercolor=[:red],xlabel="g1", ylabel="g2")
 p1 = scatter!( G[vec(IsFail_all.==0),1], G[vec(IsFail_all.==0),2],labels="safe" ,markercolor=[:blue],xlabel="g1", ylabel="g2")

 p2 = scatter( G[vec(IsFail_all),1], G[vec(IsFail_all),3],labels="fail" ,markercolor=[:red],xlabel="g1", ylabel="g3")
 p2 = scatter!( G[vec(IsFail_all.==0),1], G[vec(IsFail_all.==0),3],labels="safe" ,markercolor=[:blue],xlabel="g1", ylabel="g3")
 
 p3 = scatter( G[vec(IsFail_all),2], G[vec(IsFail_all),3],labels="fail" ,markercolor=[:red],xlabel="g2", ylabel="g3")
 p3 = scatter!( G[vec(IsFail_all.==0),2], G[vec(IsFail_all.==0),3],labels="safe" ,markercolor=[:blue],xlabel="g2", ylabel="g3")
 
 plot(p1, p2, p3, layout = l)

