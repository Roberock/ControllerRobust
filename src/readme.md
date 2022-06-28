## Functions for data analysis and to run the solver


```Julia
function g_controller(design::AbstractVector{<:Real}, 
    θ::Array,
    tvec::StepRangeLen=0:0.002:25,
    do_linear::Bool=true)

    # INPUTS 
    # design - the controller gains design = [coeffs numerator ceffs denominator]
    # θ - system parameters (uncertain parameters) 
    # tvec - discretized time
    # do_linear - [true,false]
    
    # OUTPUTS: the reliabiliy performance funtions of the controller
    # g[1] : hurwitz stability
    # g[2] : control effort 
    # g[3] : settling time 

    Ap, Bp, Cp, Dp = build_plant(θ); # multi state matrices 
    Bp1=Bp[:,1];
    Bp2=Bp[:,2]; 
    # define transfer function for the controller
    Ac, Bc, Cc, Dc = tf2ss(design[1:4], design[5:end]) 
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
        y , t =impulse(clsys,tvec);
        u, t  = lsim(ss(Ac, Bc, Cc', Dc),y,t);    
    else
        # not yet available: implement non-linear solver 
        y , t =impulse(clsys,tvec);
        u, t  = lsim(ss(Ac, Bc, Cc', Dc),y,t);    
    end
    pos_t= findfirst(t .>15)
    maxy=0.1;     
    g2=max(abs(y[pos_t]))-maxy;   # control effort 
    maxu=0.5;    
    g3=findmax(abs.(u))[1]-maxu;  # settling time 

    return g1, g2, g3 ,y
end 
```


```Julia
function EstimateReliabilityScores(G::Array) 
    ## G are the samples of the reliability functions 
    Ng=size(G,2); # number of requirements 
    Sev_ind=zeros(1,Ng); #
    W = maximum(G',dims=1);  # w=max(g1,g2,...,gNg) worst-case reliability score 
    IsF_all=W.>=0; # failure indicator for any failure
    IsF_ind= G.>=0;  # failure indicator for individual failures
    #  failure probability 
    Pf_ind = mean(G.>=0,dims=1);    # for the individual requirements
    Pf_all = mean(W.>=0);          # for the  combined  requirements
   
    for i=1:Ng
       Idxs = findall(G[:,i].>0);
       if isempty(Idxs)
        Sev_ind[i] = 0; 
       else
        Sev_ind[i] = mean(G[Idxs,i]);
       end
    end 
    Sev_all =mean(W[findall(W.>0)]);

    return Pf_ind, Pf_all, Sev_ind, Sev_all, W, IsF_all, IsF_ind
end
```

```Julia
function MonteCarloController( 
    dnom::AbstractVector{<:Real}, 
    θ::Array,
    tvec::StepRangeLen=0:0.002:25,
    do_linear::Bool=true) 
    ## function to evaluate samples of the uncertainty for a design d  
    Nsamples=size(θ,1) 
    G=zeros(Nsamples,3);
    Y=zeros(Nsamples,length(tvec));
    for i =1:Nsamples 
        g1, g2, g3, y = g_controller(dnom,θ[i,:],tvec, do_linear) 
        G[i,:]= [g1 g2 g3];
        Y[i,:]= y;

        if isinteger(i/100)
            println("Number of iterations completed $i/$Nsamples")
        end
    end     
    return G, Y
end
```
