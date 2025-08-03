
function [SOInd,total_SO]=getSO(reporpop)

reporpop_cost=[reporpop.Cost];
costs=mapminmax(reporpop_cost,0,1);

S=numel(reporpop);
SO=zeros(1,S);

Z_origin=[0;0];

for i=1:S
    SO(1,i)=sqrt(sum((costs(:,i)-Z_origin).^2)); % the smaller the better 
end

[~,SOInd]=sort(SO,'ascend'); % the smaller the better
total_SO=sum(SO);

end
