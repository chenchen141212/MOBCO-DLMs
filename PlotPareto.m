
% Project Title: Multi-Objective Particle Swarm Optimization 

function PlotPareto(pop,rep)

%% find the optimal solution point
    store=zeros(size(rep,1),1); % store the FMeasure values        
    for i=1:size(rep,1)
       store(i,1)=rep(i).Out.Acc; 
    end    
    [~,index]=max(store);     % the maimal F value is the fianl clustering result
%% plot pop and rep
%     pop_costs=[pop.Cost];
%     rep_costs=[rep.Cost];   
%     total_costs=[pop_costs rep_costs];
%     costs=mapminmax(total_costs,0,1); 
%     plot(costs(1,1:size(pop_costs,2)),costs(2,1:size(pop_costs,2)),'bo');
%     hold on;    
%     plot(costs(1,size(pop_costs,2)+1:end),costs(2,size(pop_costs,2)+1:end),'r*');
 
%% plot rep
    rep_costs=[rep.Cost];
    costs=mapminmax(rep_costs,0,1);
    plot(costs(1,:),costs(2,:),'b*');
    hold on;
    plot(costs(1,index),costs(2,index),'r*');
    hold on; 

 %% plot x and y   
    axis([-0.1 1.1, -0.1 1.1]);
    xlabel('Compactness','FontSize',15);
    ylabel('Separateness^{v}','FontSize',15)   
    hold on;    
    grid on;    
    hold off;

end