% This function is used to delete the repeated rep member

function popNodominated=DeleteSamePopRepMember(pop,rep)

    popNodominated = pop(~[pop.IsDominated]); % find the nondominated pop

    pop_Costs=[popNodominated.Cost];
    
    n1=size(rep,1);
    n2=size(pop_Costs,2);
    
    % delete the same members of rep and pop
    for i=1:n1
        mark1 =[]; 
        for j=1:n2            
            if rep(i).Cost== pop_Costs(:,j)
                mark1 = [mark1,j];
            end
        end
        popNodominated(mark1) = [];
        pop_Costs(:,mark1) = [];
        n2=size(pop_Costs,2);
    end
       
end

