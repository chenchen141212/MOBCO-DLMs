function   pop=poprank(pop,nPoP)
    
    [minf,~]=min([pop.Cost],[],2);
    
    for i=1:nPoP      
        PopCost(i)=sqrt(sum((pop(i).Cost-minf).^2)); % the smaller the better 
    end

    [~,PopInd]=sort(PopCost);
    for s=1:nPoP
        pop(PopInd(s)).ReporPopInd=s;
    end
    
end