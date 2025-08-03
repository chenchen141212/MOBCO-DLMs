% this function is used to deal with invalid clusters (data points <=3)

% refer <Automatic kernel clustering with a Multi-Elitist Particle Swarm
% Optimization Algorithm> <Elastic Differential Evolution for Automatic
% Data Clustering>

function [m, data_ind] =invalidcluster_mixed(data, m, idx_num, idx_cat)
    
    n=size(data,1);                             % number of data samples    
    K=size(m,1);                                % number of clusters
    avernumber=floor(n/K);                      % the number of data in K-1 clusters
    
    opind=zeros(n,1);                           % optimal ind solution
    nid=1:n;
    
    store_sn=[];
    temp=nid;
        
    for i=1:K-1        
        sn=nid(temp(randperm(n,avernumber)));   % selected number:generate data from n (data number is avernumber)      
        store_sn=[store_sn,sn];                 % store the selected number
        opind(sn)=i;                            % op: store the new 
        temp=setdiff(nid,store_sn);             % intersection
        n=size(temp,2);
    end
               
    opind(opind==0)=K;
    data_ind=[data opind];
    
    % compute the cluster centers
    for i=1:K
        m(i,idx_num)=mean(data(opind==i,idx_num)); % mean value
        m(i,idx_cat)=mode(data(opind==i,idx_cat)); % mode value
    end
    
end