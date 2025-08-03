
% This function is used to delete the repeated pop members
function pop=DeleteSamePopMember(pop)

    n=size(pop,1);

    for i=1:n-1
        mark = [];
        for j=(i+1):n 
            if pop(i).Cost==pop(j).Cost
                mark = [mark,j];
            end
        end
        pop(mark)=[];
        n=size(pop,1);
    end
     
end

