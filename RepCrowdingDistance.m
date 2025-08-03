
% compute the crowding distance of rep
function CDInd=RepCrowdingDistance(rep)
       
     Costs=[rep.Cost];
     nObj=size(Costs,1);     
     S=numel(rep);
     d=zeros(S,nObj);    
     CD=zeros(S,1);
          
    for j=1:nObj
        
        [cj, so]=sort(Costs(j,:));         
        d(so(1),j)=inf;
        d(so(end),j)=inf;
                
        % compute the K neighbor crowding distance
        for i=2:S-1
            d(so(i),j)=abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));            
        end
        
    end
    
    CD(:,1)=sum(d,2);
    [~,CDInd]=sort(CD,'descend');
    
end