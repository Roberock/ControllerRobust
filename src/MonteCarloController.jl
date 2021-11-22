
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