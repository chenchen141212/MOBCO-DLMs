% compute the entropy and purity 
% CP：算法聚类与实际类别得到的数据的交集，P和C的交集,C*P, Cid属于预测的，Pid属于UCI的
% Ci 算法聚类得到的每个类别的总数

function [Entropy, Purity]=EnAndPur(CP,Ci)

% get the number of rows and columns
[rn, cn]=size(CP); 

%% get the entropy
% compute precision
for i=1:rn
    for j=1:cn
        precision(i,j)=CP(i,j)/Ci(1,i);    
    end
end

% compute ei(i,j)
for i=1:rn
    for j=1:cn
        ei(i,j)=precision(i,j)*log2(precision(i,j));    
    end
end

%compute ei_sum
for i=1:rn
    ei_sum(i)=-nansum(ei(i,:)); % nansum - deal with matrix includes NaN
end

% compute mi*ei_sum(i)
for j=1:rn  %原始代码：j=1:cn
    mmi(j)=Ci(1,j)*ei_sum(j);
end

% compute entropy
Entropy=nansum(mmi)/nansum(Ci);

%% compute Purity
%找出最大的一类
for i=1:rn
     pr_max(i)=max(precision(i,:));    
end
%计算类别数量
for j=1:rn %原始代码：j=1:cn
    nni(j)=Ci(1,j)*pr_max(j);
end
Purity=nansum(nni)/nansum(Ci);
end