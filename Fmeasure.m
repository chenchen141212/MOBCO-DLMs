
function [FMeasure,Accuracy,Entropy,Purity] = Fmeasure(P,C)
% P: UCI datalabel
% C: predicted datalabel

S = length(C);       % 样本总数
pp = unique(P);      %% 返回与P中相同的数据，但是不包含重复项。pp已排序
cc = unique(C);      %% 返回与C中相同的数据，但是不包含重复项。cc已排序
P_size = length(pp); % 人工标记的簇的个数,即类别数
C_size = length(cc); % 算法计算的簇的个数，即类别数

% Pid,Cid：非零数据：第i行非零数据代表的样本属于第i个簇
Pid = double(ones(P_size,1)*P' == pp.*ones(1,S) );   %double为双精度
Cid = double(ones(C_size,1)*C' == cc*ones(1,S) );    %double为双精度
CP = Cid*Pid';      % P和C的交集,C*P, Cid属于预测的，Pid属于UCI的
Pj = sum(CP,1);     % 按照列求和，结果为行向量，P在C各个簇中的个数
Ci = sum(CP,2);     % 按照行求和，结果为列向量，C在P各个簇中的个数
 
precision = CP./( Ci*ones(1,P_size) );
recall = CP./( ones(C_size,1)*Pj );
F = 2*precision.*recall./(precision+recall);

% 得到一个总的F值
FMeasure = sum( (Pj./sum(Pj)).*max(F) );
Accuracy = sum(max(CP,[],2))/S;

%得到聚类效果 Entropy和Purity
C_i=Ci';
[e1,p1]=EnAndPur(CP ,C_i);
Entropy=e1;
Purity=p1;

end