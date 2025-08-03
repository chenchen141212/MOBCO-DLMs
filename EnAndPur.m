% compute the entropy and purity 
% CP���㷨������ʵ�����õ������ݵĽ�����P��C�Ľ���,C*P, Cid����Ԥ��ģ�Pid����UCI��
% Ci �㷨����õ���ÿ����������

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
for j=1:rn  %ԭʼ���룺j=1:cn
    mmi(j)=Ci(1,j)*ei_sum(j);
end

% compute entropy
Entropy=nansum(mmi)/nansum(Ci);

%% compute Purity
%�ҳ�����һ��
for i=1:rn
     pr_max(i)=max(precision(i,:));    
end
%�����������
for j=1:rn %ԭʼ���룺j=1:cn
    nni(j)=Ci(1,j)*pr_max(j);
end
Purity=nansum(nni)/nansum(Ci);
end