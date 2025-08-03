
function KD = KprototypeDistance(data, centroids,idx_num,idx_cat,mu)
    
    m = size(data,1);
    K = size(centroids, 1);
    
    % number of numerical attributes
    p = size(idx_num, 2);

    % kprototype distance between datasets and centroids
    KD=zeros(m,K);

    % estimate a good value for gamma if no value is provided (see Huang[1997])
    if ~isempty(idx_num) % mixed attributes
        if mu == ' ' % mu is a weight for categorical attributes in the cluster l
            mu = 0.5 * std(reshape(data(:,idx_num), 1, m*p));
        end
    else
        mu=1; % all attributes are categorical attributes
    end

    for i=1:m  % m is the number of data samples
      for  j=1:K 
            KD(i,j) = euclideanDissim(data(i,idx_num), centroids(j,idx_num)) + mu * matchingDissim(data(i,idx_cat), centroids(j,idx_cat));
      end
    end

end

% euclidean dissimilarity between a vector and the centroids matrix
function eucDissim = euclideanDissim(x, centroids)
    eucDissim = sum((x - centroids).^2, 2);
end

% simple matching dissimilarity between a vector and the centroids matrix
function matDissim = matchingDissim(x, centroids)
    matDissim = sum(x ~= centroids, 2);
end


