
function rep = Perturbation_rep(rep,CostFunction,VarMax,VarMin,idx_num,idx_cat,cn,data)
        
    repold=rep;
    
    if isempty(idx_num) % when idx_num=0
        temp_rep=idx_cat(1,randperm(size(idx_cat,2),1)); 
         for i=1:numel(rep)
            for jj=1:cn
                tempindex=randperm(size(data,1),1);
                rep(i).Position(jj,temp_rep)=data(tempindex,temp_rep);
            end
             rep(i)= boundarylimitations(rep(i),idx_num,idx_cat,VarMax,VarMin);
             [rep(i).Cost, rep(i).Out]=CostFunction(rep(i).Position);
         end     
    elseif isempty(idx_cat) % when idx_cat=0
        temp_rep=idx_num(1,randperm(size(idx_num,2),1));
         for i=1:numel(rep)
             rep(i).Position(:,temp_rep)=VarMin(1,temp_rep)+(VarMax(1,temp_rep)-VarMin(1,temp_rep)).*rand(cn,1);
             rep(i)= boundarylimitations(rep(i),idx_num,idx_cat,VarMax,VarMin);
             [rep(i).Cost, rep(i).Out]=CostFunction(rep(i).Position);
         end
    else % when idx_num~=0 and idx_cat~=0
        temp1_rep=idx_num(1,randperm(size(idx_num,2),1));
        temp2_rep=idx_cat(1,randperm(size(idx_cat,2),1)); 
         for i=1:numel(rep)
             rep(i).Position(:,temp1_rep)=VarMin(1,temp1_rep)+(VarMax(1,temp1_rep)-VarMin(1,temp1_rep)).*rand(cn,1);
             for jj=1:cn
                tempindex=randperm(size(data,1),1);
                rep(i).Position(jj,temp2_rep)=data(tempindex,temp2_rep);
             end
             rep(i)= boundarylimitations(rep(i),idx_num,idx_cat,VarMax,VarMin,data);
             [rep(i).Cost, rep(i).Out]=CostFunction(rep(i).Position);
         end
    end
    rep=[repold;rep];

%    rep=DetermineDomination([repold;rep]);
%      rep=DeleteSameRepMember(rep); 
%      rep=rep([rep.IsDominated]==0);
end