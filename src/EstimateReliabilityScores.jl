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
