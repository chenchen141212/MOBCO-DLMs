
% Pos is the position, data is the data set   
function [z, out] = ClusteringCost_mixed(Pos,data,output,idx_num, idx_cat)     

    %% compute the cost
    [Pos,IntraD,InterD_new,out]=ObjectiveFunction_mixed(data,Pos,idx_num,idx_cat); %% 0(N_data x cn x nVar)+0(cn xnVar+cn logcn)
    out.m=Pos;
    predictedind=out.ind;
    
    z1=IntraD;
    z2=InterD_new;
        
    z=[z1,z2]';
    
    %% compute AC,ARI,PRE,RE
    [Acc,ARI,Pre,Fm]=ComputeEvaluationIndex_mixed(predictedind,output,data,Pos,idx_num, idx_cat); % O(cn² × N_data)。
    out.Acc=Acc;
    out.ARI=ARI;
    out.Pre=Pre;
    out.Fm=Fm;

end