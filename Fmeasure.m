
function [FMeasure,Accuracy,Entropy,Purity] = Fmeasure(P,C)
% P: UCI datalabel
% C: predicted datalabel

S = length(C);       % ��������
pp = unique(P);      %% ������P����ͬ�����ݣ����ǲ������ظ��pp������
cc = unique(C);      %% ������C����ͬ�����ݣ����ǲ������ظ��cc������
P_size = length(pp); % �˹���ǵĴصĸ���,�������
C_size = length(cc); % �㷨����Ĵصĸ������������

% Pid,Cid���������ݣ���i�з������ݴ�����������ڵ�i����
Pid = double(ones(P_size,1)*P' == pp.*ones(1,S) );   %doubleΪ˫����
Cid = double(ones(C_size,1)*C' == cc*ones(1,S) );    %doubleΪ˫����
CP = Cid*Pid';      % P��C�Ľ���,C*P, Cid����Ԥ��ģ�Pid����UCI��
Pj = sum(CP,1);     % ��������ͣ����Ϊ��������P��C�������еĸ���
Ci = sum(CP,2);     % ��������ͣ����Ϊ��������C��P�������еĸ���
 
precision = CP./( Ci*ones(1,P_size) );
recall = CP./( ones(C_size,1)*Pj );
F = 2*precision.*recall./(precision+recall);

% �õ�һ���ܵ�Fֵ
FMeasure = sum( (Pj./sum(Pj)).*max(F) );
Accuracy = sum(max(CP,[],2))/S;

%�õ�����Ч�� Entropy��Purity
C_i=Ci';
[e1,p1]=EnAndPur(CP ,C_i);
Entropy=e1;
Purity=p1;

end