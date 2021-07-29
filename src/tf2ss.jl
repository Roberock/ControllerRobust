using LinearAlgebra

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