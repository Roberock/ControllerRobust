
function MonteCarloController(
    dnom::AbstractVector{<:Real},    θ::Array) 
    ## function to evaluate samples of the uncertainty for a design d

    Nsamples=size(θ,1)
    G=zeros(Nsamples,3);
    for i =1:Nsamples 
        do_linear=true
        g1, g2, g3 = g_controller(dnom,θ[i,:],do_linear) 
        G[i,:]= [g1 g2 g3];

        if isinteger(i/100)
            println("Number of iterations completed $i/$Nsamples")
        end
    end     
    return G
end