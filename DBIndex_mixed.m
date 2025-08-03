
% compute the DB index
function DB= DBIndex_mixed(data, m, idx_num, idx_cat)

    % Calculate Distance Matrix
    KD=KprototypeDistance(data, m, idx_num, idx_cat, ' ');
    d=KD;
    
    % Assign Clusters and Find Closest Distances
    [mindd, ind] = min(d, [], 2);
    
    nk=zeros(1,size(m,1));
    for k=1:size(m,1)
        nk(k)=sum(ind==k);         % number of data in a cluster
    end
    
    % compute the intra cluster distance
    temp1=zeros(1,size(m,1));
    for kk=1:size(m,1)
        temp1(1,kk)=sum(mindd(ind==kk).^2)/nk(kk);
    end

    % compute the inter cluster distance
    MMD=KprototypeDistance(m,m,idx_num, idx_cat, ' ');
    
    r = zeros(size(m,1));
    for i=1:size(m,1)-1
        for j=i+1:size(m,1)
            r(i,j) = (temp1(i)+temp1(j))/MMD(i,j);
            r(j,i) = r(i,j);
        end
    end
  
    R=max(r);

    DB = mean(R);
       
end