
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