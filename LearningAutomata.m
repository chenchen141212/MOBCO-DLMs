function cate_infor=LearningAutomata(individual,idx_cat,newcost,lastcost,cate_infor)
    
    alpha=0.1;
    beta=0.01;
    
    if newcost <= lastcost
        for j=1:size(idx_cat,2) % the jth categorical attributes for a bacteriaum
            temp=individual(:,idx_cat(j));
            temp=unique(temp);
            row=zeros(numel(temp),1);
            for n=1:numel(temp)
                [row(n),~]=find(temp(n)==cate_infor(j).avaliable_total);
                cate_infor(j).Prob2(row(n))=cate_infor(j).Prob2(row(n))+alpha.*(1 - cate_infor(j).Prob2(row(n)));  % desired action
            end
            cate_infor(j).Prob2(cate_infor(j).nj~=row)=(1-alpha)*cate_infor(j).Prob2(cate_infor(j).nj~=row);
        end
    else
        for j=1:size(idx_cat,2) % the jth categorical attributes for a bacteriaum
            temp=individual(:,idx_cat(j));
            temp=unique(temp);
            row=zeros(numel(temp),1);
            for n=1:numel(temp) 
                [row(n),~]=find(temp(n)==cate_infor(j).avaliable_total);
                cate_infor(j).Prob2(row(n))=(1-beta).*cate_infor(j).Prob2(row(n));  %  non-desired action
            end
            cate_infor(j).Prob2(cate_infor(j).nj~=row)=(beta/(cate_infor(j).nj-1))+(1-beta)*cate_infor(j).Prob2(cate_infor(j).nj~=row);
        end                   
    end
    
end