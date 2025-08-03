%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com

%%
% pop means rep here, nondominated solutions structure
% nGrid=7: Number of Grids per Dimension
% alpha=0.1: Inflation Rate

function Grid=CreateGrid(pop,nGrid,alpha)

    c=[pop.Cost];

    cmin=min(c,[],2);    % row minmal value
    cmax=max(c,[],2);
    
    dc=cmax-cmin;
    cmin=cmin-alpha*dc;
    cmax=cmax+alpha*dc;
    
    nObj=size(c,1);
    
    empty_grid.LB=[];
    empty_grid.UB=[];
    Grid=repmat(empty_grid,nObj,1);
    
    for j=1:nObj
        
        cj=linspace(cmin(j),cmax(j),nGrid+1); % linspace: generate lineraly spaced vector
        
        Grid(j).LB=[-inf cj];
        Grid(j).UB=[cj +inf];
        
    end

end