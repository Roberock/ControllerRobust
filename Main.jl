#using Pkg
# Pkg.add("DifferentialEquations"); Pkg.add("ControlSystems"); kg.add("Plots")
# import Pkg; Pkg.add("ControlToolbox") 
using ControlSystems, StatsPlots
using Distributed, StatsBase, Distributions 

include(string(pwd(), "\\src\\build_plant.jl"))
include(string(pwd(), "\\src\\g_controller.jl"))
include(string(pwd(), "\\src\\tf2ss.jl"))
include(string(pwd(), "\\src\\MonteCarloController.jl"))

 
##### START building the controller

## % Parameter θ
θref=[1, 1, 1, 0, 1,0.001]      ## nominal θ
names=["m1" "m2" "kl" "kn" "λ" "τ" ]     ## names of parameters
# admissible ranges
θlims=[0.1    2;        0.1    2;
       0.05   1.75;     -1     1;
       0.2    1.8;      0.0001 0.3];             
       
# nominal design      
dnom=[-0.1324, 0.3533, 0.6005, 0.0728, 0.5503, 1.4175, 2.6531, 2.4802, 1];  

### build the system state-space matrix 
A, B, C, D = build_plant(θref); 

### define converter from trasfer function to system  
Ac, Bc, Cc, Dc=tf2ss(dnom[1:4],dnom[5:end]);  
###  define controller reliability function 
 
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

## load samples
using CSV, DataFrames

File_path=string(pwd(), "\\Theta_small.csv");
df = CSV.read(File_path,DataFrame);
Theta_samples=Matrix(df);

G_samples  =  MonteCarloController(dnom,Theta_samples)
Pf, Severity = EstimateReliabilityScores(G_samples)
corrplot(G_samples,labels=["g$i" for i=1:3])

 

