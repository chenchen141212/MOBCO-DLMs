
% This function is used to delete the repeated rep members
function rep=DeleteSameRepMember(rep)

    n=size(rep,1);

    for i=1:n-1
        mark = [];
        for j=(i+1):n 
            if rep(i).Cost==rep(j).Cost
                mark = [mark,j];
            end
        end
        rep(mark)=[];
        n=size(rep,1);
    end
     
end

