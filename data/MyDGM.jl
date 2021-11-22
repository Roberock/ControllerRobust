function MyDGM(θref, θlims::Array, Nsamples) 
    ## function to evaluate samples of the uncertainty for a design d  
     
    
    θmyDGM=  θref.*(rand(6,Nsamples).-0.5);
    for i =1:size(θlims,1) 
        θmyDGM[i,:]=min.(θmyDGM[i,:],θlims[i,2]);
        θmyDGM[i,:]=max.(θmyDGM[i,:],θlims[i,1]);
    end

    return θmyDGM
end