
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)

function pop=DetermineDomination(pop)

    nPop=numel(pop);

    % initialization:let every bacteria is non-dominated first
    for i=1:nPop
        pop(i).IsDominated=false; 
    end

    % determine the dominated or non-dominated
    for i=1:nPop-1
        for j=i+1:nPop
            if Dominates(pop(i),pop(j))  % pop(i) dominated pop(j)
               pop(j).IsDominated=true;
            end
            
            if Dominates(pop(j),pop(i))  % pop(j) dominated pop(i)
               pop(i).IsDominated=true;
            end
            
        end
    end  
end