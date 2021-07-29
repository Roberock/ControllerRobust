#using Pkg
# Pkg.add("DifferentialEquations"); Pkg.add("ControlSystems"); kg.add("Plots")
# import Pkg; Pkg.add("ControlToolbox") 
using ControlSystems 
using LinearAlgebra, PyPlot   
##### START building the controller

## % Parameter θ
θ=[1 1 1 0 1 0.001]      ## nominal θ
names=["m1" "m2" "kl" "kn" "λ" "τ" ]     ## names of parameters
# admissible ranges
θlims=[0.1    2;        0.1    2;
       0.05   1.75;     -1     1;
       0.2    1.8;      0.0001 0.3];             
       
# nominal design      
dnom=[-0.1324 0.3533 0.6005 0.0728 0.5503 1.4175 2.6531 2.4802 1];  

###
function  build_plant(θ::AbstractVector{<:Real})
    p=θ
    # linearization about zero (not about stable oscillator). 
    # Nonlinear spring does not show up
        if length(θ)==2 
            p_out=[p[1] 1 p[2] 0 1 0.001];
        else  
            p_out=p;
           #  p_out[1]=p[1];    # m1  
           # p_out[2]=p[2];    # m2
           # p_out[3]=p[3];    # klinear
           # p_out[4]=p[4];    # knonlin
           # p_out[5]=p[5];    # control effectivenes
           # p_out[6]=p[6];    # time delay
        end
    ratio=1/p_out[6];
    A=[0 0 1 0 0;
             0 0 0 1 0;
             -p_out[3]/p_out[1]    p_out[3]/p_out[1]   0   0   p_out[5]/p_out[1]; 
             p_out[3]/p_out[2]     -p_out[3]/p_out[2]  0   0     0; 
              0                   0                    0   0    -ratio];
              
    Bp1=[0 0 0 0 ratio]';
    Bp2=[0 0 0 1/p_out[2] 0]';
    B=[Bp1 Bp2];
    C=[0 1 0 0 0];
    D=[0 0];
    x0=Bp2;     #simulating impulse: zero input and IC=Base

    return A, B, C, D, x0, p_out
end
 
A, B, C, D = build_plant(θ) 
### define converter from trasfer function to system

 
 
function tf2ss(Num::AbstractVector{<:Real},Den::AbstractVector{<:Real})
    # TF2SS  Transfer function to state-space conversion.
    #   calculates the state-space representation: 
    #       x = Ax + Bu
    #       y = Cx + Du 
    #   of the system:
    #               NUM(s)
    #       H(s) = --------
    #               DEN(s) 
    if isempty(Num) && isempty(Den)
        a = [];        b = [];
        c = [];        d = [];
    else 
        denStrip=Den; 
        nnum= length(Num); 
        nden = length(Den);  
        # Check for proper numerator
        if (nnum > nden) 
            # Try to strip leading zeros to make proper
            numStrip = Num[(nnum-nden+1):nnum];
        else 
            numStrip = [ zeros(1,nden-nnum); Num];
        end  
        numNormalized = numStrip/denStrip[1];
        denNormalized = denStrip/denStrip[1]; 
        d = numNormalized[1,:];
        c = numNormalized[2:nden] - numNormalized[1]*denNormalized[2:nden];
         

        if nden == 1
            a = [];     b = [];            c = [];
        else
            a = vcat(-denNormalized[2:nden]' , Matrix{Float64}(I, nden-2, nden-1));
            b = Matrix{Float64}(I, nden-1, 1) ; 
        end
    end 
    return a,b,c,d
end
 
Ac, Bc, Cc, Dc=tf2ss(dnom[1:4],dnom[5:end])  
 ##  define controller reliability function 

function g_controller(design,θ,do_linear::Bool=true)

    # INPUTS 
    # design - the controller gains design = [coeffs numerator ceffs denominator]
    # θ - system parameters (uncertain parameters) 
    # do_linear - [true,false]
    
    # OUTPUTS: the reliabiliy performance funtions of the controller
    # g[1] : hurwitz stability
    # g[2] : control effort 
    # g[3] : settling time 
    Ap, Bp, Cp, Dp = build_plant(θ); # multi state matrices 
    Bp1=Bp[:,1];
    Bp2=Bp[:,2]; 
    # define transfer function for the controller
    Ac, Bc, Cc, Dc = tf2ss(design[1:4],design[5:end]) 
    com_poles=eigvals(Ac); 
     # define system state-space model
    Acl=vcat( hcat(Ap-Bp1.*Dc*Cp, Bp1*Cc'), hcat( -Bc*Cp , Ac));
    Bcl=vcat( Bp2 , zeros(size(Ac,1),1));
    Ccl=hcat(Cp,   zeros(1,size(Ac,1)));
    Dcl=0;

    clsys=ss(Acl,Bcl,Ccl,Dcl); 
    sys_poles=pole(clsys) # system poles
    numclpoles= length(sys_poles)  #  number of poles
    poles=[sys_poles;com_poles]; # poles of closed loop and of controller
    g1=max(real(poles)...);  # compute hurwitz stability
  
    
    if do_linear 
        y , t =impulse(clsys,0:0.001:25);
        u, t  = lsim(clsys,y,t);    
    else
        # to do: implement non linear solver 
        y , t =impulse(clsys,0:0.001:25);
        u, t  = lsim(clsys,y,t);    
    end
    pos_t= findfirst(t .>15)
    maxy=0.1;     g2=max(abs(y[pos_t]))-maxy;
    maxu=0.5;     g3=findmax(abs.(u))[1]-maxu;

    return g1, g2, g3
end


gnominal=g_controller(dnom,θ,true)
## Generate random samples within the limits of θ

function MonteCarloController(
    dnom::AbstractVector{<:Real},
    θref::AbstractVector{<:Real},
    Nsamples::Integer=100)

    Gsave1=[];     Gsave2=[];    Gsave3=[];
    for i =1:Nsamples 
        do_linear=true
        g1, g2, g3 = g_controller(dnom,abs.(θref.+randn(1,6)/1000),do_linear) 
        append!(Gsave1,g1 )  
        push!(Gsave2,g2 )  
        push!(Gsave3,g3 )
    end     
    return Gsave
end

Nsamples = 100;
Gsave  =  MonteCarloController(dnom,θ,Nsamples)