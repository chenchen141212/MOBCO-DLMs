
function [m,IntraD,InterD_new,out]=ObjectiveFunction_mixed(data, m, idx_num, idx_cat)

    % Calculate Distance Matrix
    KD=KprototypeDistance(data, m, idx_num, idx_cat, ' '); %(O(N_data×cn×nVar)​)
   
    % Assign Clusters and Find Closest Distances
    [mindd, ind] = min(KD, [], 2);   % ind is the cluster labels
    
    data_ind=[data ind];
    
    nk=zeros(1,size(m,1));
    
    for k=1:size(m,1)
        nk(k)=sum(ind==k);         % number of data in a cluster
    end

    if any(nk<=3)  %-------------------------------------------------------
        [m, data_ind] =invalidcluster_mixed(data, m, idx_num, idx_cat);  %% O(N_data*nVar)​ 
        ind_=data_ind(:,end);
        
        % compute the intra cluster distance(minimize)
        temp1=zeros(1,size(m,1));
        temp3=zeros(1,size(m,1));
        nkk=zeros(1,size(m,1));
        
        for kk=1:size(m,1)
            nkk(kk)=sum(ind_==kk); 
            KD_=KprototypeDistance(data_ind(ind_==kk,1:end-1), m(kk,:),idx_num, idx_cat, ' ');
            temp1(1,kk)=sum(KD_)/nkk(kk);    % related to IntraD
        end
        
        % compute the inter cluster distance(maxmize)
        temp2=zeros(1,size(m,1));            % related to InterD
        MMD=KprototypeDistance(m,m,idx_num, idx_cat, ' ');
        sortedMMD = sort(MMD,2);             % ind is the cluster labels
        
        for kk=1:size(m,1)
            temp2(1,kk)=sortedMMD(kk,2);
        end
        
    else %-----------------------------------------------------------------
        % compute the intra cluster distance
        temp1=zeros(1,size(m,1)); % related to IntraD
        for kk=1:size(m,1)
            temp1(1,kk)=sum(mindd(ind==kk))/nk(kk);
        end
        
        % compute the inter cluster distance
        temp2=zeros(1,size(m,1)); % related to InterD
        MMD=KprototypeDistance(m,m,idx_num, idx_cat, ' ');
        sortedMMD = sort(MMD,2);  % ind is the cluster labels
        for kk=1:size(m,1)
            temp2(1,kk)=sortedMMD(kk,2);
        end 
                        
    end  %-----------------------------------------------------------------
       
    %% output
    IntraD=sum(temp1); 
    InterD=sum(temp2);
    InterD_new=1/InterD;    
    
    out.m=m;
	out.ind=data_ind(:,end);
    out.IntraD=IntraD;
    out.InterD_new=InterD_new;
end