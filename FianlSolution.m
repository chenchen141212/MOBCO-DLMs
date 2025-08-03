% this function is used to select the final clustering results

function FianlOut=FianlSolution(rep)

    store=zeros(size(rep,1),1); % store the FMeasure values
        
    for i=1:size(rep,1)
       store(i,1)=rep(i).Out.Acc; 
    end
    
    [acmax,index]=max(store);     % the maimal F value is the fianl clustering result  
    FianlOut.Acc=acmax;                  % max
    FianlOut.ARI=rep(index).Out.ARI;
    FianlOut.Pre=rep(index).Out.Pre;
    FianlOut.Fm=rep(index).Out.Fm;
end