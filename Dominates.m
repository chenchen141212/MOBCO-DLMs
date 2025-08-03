
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)

function b=Dominates(x,y)

    if isstruct(x)
        x=x.Cost;
    end
    
    if isstruct(y)
        y=y.Cost;
    end

    b=all(x<=y) && any(x<y);

end