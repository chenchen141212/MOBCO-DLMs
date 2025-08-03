
function rep=getRep(pop,rep,nRep,nGrid,alpha)

    % update non-dominated repository----------------------------------
    pop=DetermineDomination(pop);  % select non-dominated from pop
    
    % delete same members in pop
    pop=DeleteSamePopMember(pop);
    
    % delete same members in rep
    rep=DeleteSameRepMember(rep);
    
    % delete same members in pop and rep
    popNodominated=DeleteSamePopRepMember(pop,rep);
    
%   Add Non-Dominated Particles to REPOSITORY
    rep=[rep;popNodominated];
    
    % Determine Domination of New Resository Members
    rep=DetermineDomination(rep); % select non-dominated from rep
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated]);
    
    % Update Grid
    Grid=CreateGrid(rep,nGrid,alpha);
    
    % Update Grid Indices
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteCrowdingOne(rep);
        end
    end
    
end