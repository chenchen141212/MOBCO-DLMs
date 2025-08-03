
function individual = boundarylimitations(individual,idx_num,idx_cat,VarMax,VarMin,data)
    
    K=size(individual.Position,1); % K is the cluster number

    for i=1:K
        for j=1:size(idx_num,2)
            if individual.Position(i,idx_num(j)) < VarMin(idx_num(j))|| individual.Position(i,idx_num(j)) > VarMax(idx_num(j))
                individual.Position(i,idx_num(j))=VarMin(idx_num(j))+(VarMax(idx_num(j))-VarMin(idx_num(j))).*rand;
            end
        end
        
        for j=1:size(idx_cat,2)
            if individual.Position(i,idx_cat(j)) < VarMin(1,idx_cat(j))|| individual.Position(i,idx_cat(j)) > VarMax(1,idx_cat(j))
                tempindex=randperm(size(data,1),1);
                individual.Position(i,idx_cat(j))=data(tempindex,idx_cat(j));
            end
        end  
    end
end