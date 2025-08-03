
function [Acc,ARI,Pre,Fm]=ComputeEvaluationIndex_mixed(cluster,output,data,m,idx_num, idx_cat)  

    k=numel(unique(cluster));
    Re=[cluster,output];  % cluster are the predicted lables,Ak are the ture lables
    n=size(Re,1);
    A1=Re(:,1)'; % predicted lables
    A2=Re(:,2)'; % ture lables
    cp=zeros(k,k);
    
    for i=1:k
        p1=find(A1==i);
        for j=1:k
          c1=find(A2==j);
          cp(i,j)=length(intersect(p1,c1));
        end
    end
    
    cp(k+1,:)=sum(cp);
    cp(:,k+1)=sum(cp,2);
    %------------------------------------------------------------------
    % Acc--accuracy
    a=zeros(1,k);
    for i=1:k  
        a(i)=max(cp(i,1:k));
    end
    Acc=sum(a)/n;
    %-------------------------------------------------------------------
    % ARI-adjusted rand index
    cp1=cp(1:k,1:k);
    r0=sum(sum((cp1.*(cp1-1))./2));
    cp2=cp(1:k,1+k)';
    r1=sum((cp2.*(cp2-1))./2);
    cp3=cp(1+k,1:k);
    r2=sum((cp3.*(cp3-1))./2);
    r3=(2*r1*r2)/(n*(n-1));
    ARI=(r0-r3)/(0.5*(r1+r2)-r3);

    %----------------------------------------------------------------
    % PRE--precision & RE--recall
    Pre=(sum(max(cp(1:k,1:k),[],2)./cp(1:k,k+1)))/k;
%     Re=(sum(max(cp(1:k,1:k))./cp(k+1,1:k)))/k;
    
    %----------------------------------------------------------------
    % fmeasure and entropy
    [Fm,~,~,~] = Fmeasure(cluster,output);
    
    %----------------------------------------------------------------
    % DB index
%     DB= DBIndex_mixed(data, m, idx_num, idx_cat);
    
end