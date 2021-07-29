###
function  build_plant(θ::AbstractVector{<:Real})
    p=θ
    # linearization about zero (not about stable oscillator). 
    # Nonlinear spring does not show up
        if length(θ)==2 
            p_out=[p[1], 1, p[2], 0, 1, 0.001];
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