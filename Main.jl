#using Pkg
# Pkg.add("DifferentialEquations"); Pkg.add("ControlSystems"); kg.add("Plots")
# import Pkg; Pkg.add("ControlToolbox") 
using ControlSystems, StatsPlots
using Distributed, StatsBase, Distributions 

include(string(pwd(), "\\src\\build_plant.jl"))
include(string(pwd(), "\\src\\g_controller.jl"))
include(string(pwd(), "\\src\\tf2ss.jl"))
include(string(pwd(), "\\src\\MonteCarloController.jl"))

 
##### #####  Define contorl problem  ##### #####

## Parameter θ (uncertain parameters)
θref=[1, 1, 1, 0, 1,0.001]      ## nominal θ
θnames=["mass1" "mass2" "kl" "kn" "λ" "τ" ]     ## names of the uncertain factors
# admissible ranges
θlims=[0.1    2;        0.1    2;
       0.05   1.75;     -1     1;
       0.2    1.8;      0.0001 0.3];             
       
# nominal design      
dnom=[-0.1324, 0.3533, 0.6005, 0.0728, 0.5503, 1.4175, 2.6531, 2.4802, 1];
# optimized design      
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


Design_Selection=4
if Design_Selection==1 # nominal
    dnom=[-0.1324, 0.3533, 0.6005, 0.0728, 0.5503, 1.4175, 2.6531, 2.4802, 1];  
elseif Design_Selection==2 # from min { \hat{F}^{-1}(1) }
    dnom=[ -0.103539257684432,	0.452027845175231,	0.553521292381580,	0.0731453358298691,	0.285291345952459, 0.703095280659925,	3.45687805905316,	2.49234017067302,	1.16716573966267]; # 0.717520994282477
elseif Design_Selection==3 # from min { \hat{F}^{-1}(0.95) }
    dnom=[-0.0276450664037319, 0.613233915821064, 0.528302902131644, 0.0732150560642133, 0.585865170685746, 1.35042413219186, 2.51075569540054, 2.17663965704309, 0.827018013978000];
elseif Design_Selection==4 # from  min { Pf }
    dnom=[-0.0831455337335365,	0.635784535253079,	0.775222237017883,	0.0981477664191349,	0.760033487719863,	1.94907797910329,	3.04966829435795,	2.73439891631291,	1.05935063251108]; #  gamma: 0.2218
end


d_names =["m1" "m2" "kl" "kn" "λ" "τ" ]     ## names of design parameters 
### Example: how to build the system state-space matrix?
A, B, C, D = build_plant(θref); 

### Example: how to define converter from trasfer function to system?  
Ac, Bc, Cc, Dc=tf2ss(dnom[1:4],dnom[5:end]);  

###  Define controller reliability function 
gnominal=g_controller(dnom,θref,true);

## Generate random samples within the limits of θ  
Gref  =  MonteCarloController(dnom,Array(θref'))

function EstimateReliabilityScores(
    G::Array) 
    Pf = mean(G.>0,dims=1);   
    Ng=size(G,2);
    Severity=zeros(1,Ng);
    for i=1:Ng
       Idxs = findall(G[:,i].>=0);
       if isempty(Idxs)
        Severity[i] = 0;
       else
        Severity[i] = mean(G[Idxs,i]);
       end
    end 
    return Pf, Severity
end

Pf, Severity = EstimateReliabilityScores(Gref);

## load data (samples of the uncertain factors)
using CSV, DataFrames
DataSetType =1
Nuber_of_samples=10^3
if DataSetType==1
    File_path=string(pwd(), "\\data\\Theta_small.csv"); ## small uncertainty
    df = CSV.read(File_path,DataFrame);
    Theta_samples=Matrix(df);
elseif DataSetType==2
    File_path=string(pwd(), "\\data\\Theta_large.csv"); ## large uncertainty
    df = CSV.read(File_path,DataFrame);
    Theta_samples=Matrix(df);
elseif DataSetType==3
    include(string(pwd(), "\\data\\MyDGM.jl"))## my Data-Generating-MechanismDGM
    Theta_samples=MyDGM(Nuber_of_samples);
end

## Monte carlo simulation of the nominal design 
G_samples  =  MonteCarloController(dnom,Theta_samples)

Pf, Severity = EstimateReliabilityScores(G_samples) ## Failure probability and severity  
#corrplot(G_samples,labels=["g$i" for i=1:3])

## Analyze failures and plot results  
WorsCaseg = maximum(G_samples,dims=2);  
IsFail_all=WorsCaseg.>0;
IsFail_individualRequirements= G_samples.>0;
 
 l = @layout [a ; b c]
 p1 = scatter( G_samples[vec(IsFail_all),1], G_samples[vec(IsFail_all),2],labels="fail" ,markercolor=[:red],xlabel="g1", ylabel="g2")
 p1 = scatter!( G_samples[vec(IsFail_all.==0),1], G_samples[vec(IsFail_all.==0),2],labels="safe" ,markercolor=[:blue],xlabel="g1", ylabel="g2")

 p2 = scatter( G_samples[vec(IsFail_all),1], G_samples[vec(IsFail_all),3],labels="fail" ,markercolor=[:red],xlabel="g1", ylabel="g3")
 p2 = scatter!( G_samples[vec(IsFail_all.==0),1], G_samples[vec(IsFail_all.==0),3],labels="safe" ,markercolor=[:blue],xlabel="g1", ylabel="g3")
 
 p3 = scatter( G_samples[vec(IsFail_all),2], G_samples[vec(IsFail_all),3],labels="fail" ,markercolor=[:red],xlabel="g2", ylabel="g3")
 p3 = scatter!( G_samples[vec(IsFail_all.==0),2], G_samples[vec(IsFail_all.==0),3],labels="safe" ,markercolor=[:blue],xlabel="g2", ylabel="g3")
 
 plot(p1, p2, p3, layout = l)

  