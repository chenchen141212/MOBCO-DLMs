%% implementing BFO/MOBFO for mixed data clustering alogorithms
clear;
close all;
clc;
clear;

for TP=1

    [dataname,data,output,inputType,idx_num,new_idx_num,idx_cat,cn] = TestingDataSets(TP);
    
    TotalRun=2;
    MaxIt=2;
        
    %% Problem Definition
    CostFunction=@(Pos) ClusteringCost_mixed(Pos,data,output,idx_num,idx_cat);     % Cost Function
    
    nVar=size(data,2);
    
    VarMin= min(data);    % Lower Bound of Variables
    VarMax= max(data);    % Upper Bound of Variables
    
    %% MOBFO Parameters
    nPoP=100;            % the whole Population Size
    nSub=2;              % the number of subpopulations
    SubPoP(1:nSub)=nPoP/nSub; % the population size of a subpopulation
    
    Fre=5;
    Fmig=2;
    Ns=4;
    Sr=ceil(nPoP/2);
    Pmig=0.25;
    
    c1=1;   % Personal Learning Coefficient
    c2=3;   % Global Learning Coefficient    
    
    Cmax=1.5;
    Cmin=0.1;
    cc=2;
    
    nRep=100;           % Repository Size, store the nondominated bacterium/particles
    nGrid=10;           % Number of Grids per Dimension
    alpha=0.1;          % Inflation Rate
    
    Acc_all_runs = cell(TotalRun, 1);
    ARI_all_runs = cell(TotalRun, 1);
    Pre_all_runs = cell(TotalRun, 1);
    Fm_all_runs = cell(TotalRun, 1);
    rep_all_runs= cell(TotalRun, 1);

    for run=1:TotalRun % parfor

        cate_infor=struct;
        %% Initialization
        empty_bacterial = struct;
        empty_bacterial.Position=[];
        empty_bacterial.Cost=[];
        empty_bacterial.Tumble=[];
        empty_bacterial.Delta=[];
        empty_bacterial.StepSize=[];
        empty_bacterial.Out=[];
        
        empty_bacterial.Best.Position=[];
        empty_bacterial.Best.Cost=[];
        empty_bacterial.Best.Out=[];
        
        empty_bacterial.IsDominated=[];
        empty_bacterial.GridIndex=[];
        empty_bacterial.GridSubIndex=[];
        
        empty_bacterial.ReporPopInd=[];
        
        pop=repmat(empty_bacterial,nPoP,1);
        
        for i=1:nPoP
            for j=1:size(idx_num,2)
                pop(i).Position(:,idx_num(j))=VarMin(1,idx_num(j))+(VarMax(1,idx_num(j))-VarMin(1,idx_num(j))).*rand(cn,1);
            end
            for j=1:size(idx_cat,2)        
                for jj=1:cn
                    tempindex=randperm(size(data,1),1);
                    pop(i).Position(jj,idx_cat(j))=data(tempindex,idx_cat(j));
                end
            end
            
            [pop(i).Cost, pop(i).Out]=CostFunction(pop(i).Position);
            
            % Update Personal Best
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
            pop(i).Best.Out=pop(i).Out;    
        end

        % Determine Domination
        pop=DetermineDomination(pop);
        rep=pop(~[pop.IsDominated]);
        rep=DeleteSameRepMember(rep);
        Grid=CreateGrid(rep,nGrid,alpha);
        
        for i=1:numel(rep)
            rep(i)=FindGridIndex(rep(i),Grid);
        end
        
        % obtain the total avaliable values of data
        for j=1:numel(idx_cat)
            cate_infor(j).avaliable_total =unique(data(:,idx_cat(j))); % the total avaliable values
            cate_infor(j).nj=size(cate_infor(j).avaliable_total,1);    % nj is the total number of avaliable values of jth categorical attribute                     
            cate_infor(j).Prob1=zeros(1,cate_infor(j).nj);              % the initial probability of jth categorical attribute
            cate_infor(j).Prob2=repmat(1/cate_infor(j).nj,1,cate_infor(j).nj);
            cate_infor(j).Prob=zeros(1,cate_infor(j).nj); 
        end
         
        %% MOBCO_DLMs Main Loop
        count_rep=0; % record the times that the rep didn't change

        for it=1:MaxIt
        
            % compute the total rank of every nondominated solution----------------
            RepCDInd=RepCrowdingDistance(rep);  % crowding distance
            
            [RepSOInd,total_SO]=getSO(rep);  % origin-nondominated solution distance
            
            if it==1
                total_SOH=total_SO;
            end
            
            TempRepInd=RepCDInd'+RepSOInd;
            
            [~,RepInd]=sort(TempRepInd); % ascend
            
            for s=1:numel(rep)
                rep(s).ReporPopInd=RepInd(s);
            end
        
            % compute the total rank of every population---------------------------
            pop=poprank(pop,nPoP);
            
            % build the categorical attribute candidate pool 
            temppop = repmat(pop(1), length(pop), 1); % Pre-allocate the size of the temppop
            count_pool = 0; % record the number of entities actually added
            
            for i = 1:length(pop)
                if pop(i).ReporPopInd <= SubPoP
                    count_pool = count_pool + 1; % update the counter
                    temppop(count_pool) = pop(i); % save individuals that meet the requirements
                end
            end
            
            temppop(count_pool+1:end) = [];  % delete the bottom half of the individual
        
            templearned_candidate=[rep;temppop];  
            learnedsize=size(templearned_candidate,1);
            learned_candidate=templearned_candidate(randperm(learnedsize,SubPoP(1)));
            
            % compute the occurrence probability
            for j=1:size(idx_cat,2) % the jth categorical attributes for a bacteriaum
                cate_infor(j).Count=zeros(1,cate_infor(j).nj);
                for n=1:cate_infor(j).nj
                    for s=1:SubPoP
                        for k=1:cn
                            if learned_candidate(s).Best.Position(k,idx_cat(j))== cate_infor(j).avaliable_total(n) 
                                cate_infor(j).Count(n)=cate_infor(j).Count(n)+1;
                            end
                        end
                    end                    
                end
                for n=1:cate_infor(j).nj
                    cate_infor(j).Prob1(n)=cate_infor(j).Count(n)/sum(cate_infor(j).Count(:));            
                    cate_infor(j).Prob(n)=cate_infor(j).Prob1(n)+cate_infor(j).Prob2(n);
                end
                        
                [~,betterProb]= sort(cate_infor(j).Prob,'descend');
            end

            % chemotaxis-----------------------------------------------------------
            %% for the first sub-colony (ring topology)
            for i1=1:SubPoP(1)
                    
                leader= rep([rep.ReporPopInd]==randperm(min(numel(rep),3),1)); 
                lastcost=pop(i1).Cost;   
                
                pop(i1).Delta=(2*round(rand(1,nVar))-1).*rand(1,nVar); % random walk
                pop(i1).Tumble=pop(i1).Delta/sqrt(pop(i1).Delta*(pop(i1).Delta)');        
                pop(i1).StepSize=Cmin+exp(-cc*(it/MaxIt)^2)*(Cmax-Cmin);
                
                % find neighbours
                if i1==1
                    friend1ind=SubPoP(1); 
                else
                    friend1ind=i1-1;           
                end
                
                if i1==SubPoP(1) 
                    friend2ind=1;
                else
                    friend2ind=i1+1;
                end
                
                tempmbest=struct;
                tempmbest(1).Position=pop(i1).Best.Position;
                tempmbest(2).Position=pop(friend1ind).Best.Position;
                tempmbest(3).Position=pop(friend2ind).Best.Position;
                
                tempmbest(1).Cost=pop(i1).Best.Cost;
                tempmbest(2).Cost=pop(friend1ind).Best.Cost;
                tempmbest(3).Cost=pop(friend2ind).Best.Cost;
        
                tempmbest=DetermineDomination(tempmbest);
                tempzero=find([tempmbest.IsDominated]==0);
                nbest=tempmbest(tempzero(1));
        
                % update numeric part
                if ~isempty(idx_num)
                    pop(i1).Position(:,idx_num)=pop(i1).Position(:,idx_num)+pop(i1).StepSize*(...
                        c1*rand.*(nbest.Position(:,idx_num)-pop(i1).Position(:,idx_num))+...
                        c2*rand.*(leader.Position(:,idx_num)-pop(i1).Position(:,idx_num))+pop(i1).Tumble(:,idx_num));
                end
                
                % update categorical part
                if ~isempty(idx_cat)
                    pop(i1).Position(randperm(cn,1),idx_cat(j))=cate_infor(j).avaliable_total(betterProb(1));
                end
                                
                % position limitations        
                pop(i1) = boundarylimitations(pop(i1),idx_num,idx_cat,VarMax,VarMin);
                
                [pop(i1).Cost, pop(i1).Out]=CostFunction(pop(i1).Position);        
                newcost=pop(i1).Cost;
                
                % Learning Automata strategy to categorical attributes
                if ~isempty(idx_cat)
                    cate_infor=LearningAutomata(pop(i1).Position,idx_cat,newcost,lastcost,cate_infor);
                end
                
                % Update Personal Best
                if Dominates(pop(i1),pop(i1).Best)
                    pop(i1).Best.Position=pop(i1).Position;
                    pop(i1).Best.Cost=pop(i1).Cost;
                    pop(i1).Best.Out=pop(i1).Out; 
                    
                elseif Dominates(pop(i1).Best,pop(i1))
                    % Do Nothing
                else
                    if rand<0.5
                        pop(i1).Best.Position=pop(i1).Position;
                        pop(i1).Best.Cost=pop(i1).Cost;
                        pop(i1).Best.Out=pop(i1).Out; 
                    end
                end      
                        
                % swim-------------------------------------------------------------
                m=0;
                while m<Ns
                    
                    if newcost<= lastcost              % compare cost to find better cost
                        lastcost=newcost;

                        tempmbest1=struct;
                        tempmbest1(1).Position=pop(i1).Best.Position;
                        tempmbest1(2).Position=pop(friend1ind).Best.Position;
                        tempmbest1(3).Position=pop(friend2ind).Best.Position;
        
                        tempmbest1(1).Cost=pop(i1).Best.Cost;
                        tempmbest1(2).Cost=pop(friend1ind).Best.Cost;
                        tempmbest1(3).Cost=pop(friend2ind).Best.Cost;
                        
                        tempmbest1=DetermineDomination(tempmbest1);
                        tempzero=find([tempmbest1.IsDominated]==0);
                        nbest=tempmbest1(tempzero(1));                
                        
                        % % update numeric part
                        if ~isempty(idx_num)
                            pop(i1).Position(:,idx_num)=pop(i1).Position(:,idx_num)+pop(i1).StepSize*(...
                                c1*rand.*(nbest.Position(:,idx_num)-pop(i1).Position(:,idx_num))+...
                                c2*rand.*(leader.Position(:,idx_num)-pop(i1).Position(:,idx_num)));
                        end
                        
                        % update categorical part
                        if ~isempty(idx_cat)
                           pop(i1).Position(randperm(cn,1),idx_cat(j))=cate_infor(j).avaliable_total(betterProb(1));
                        end
                        
                        % position limitations
                        pop(i1) = boundarylimitations(pop(i1),idx_num,idx_cat,VarMax,VarMin);
                        
                        [pop(i1).Cost, pop(i1).Out]=CostFunction(pop(i1).Position);    % calculate cost
                        newcost=pop(i1).Cost;
                        
                        % Learning Automata strategy to categorical attributes
                        if ~isempty(idx_cat)
                            cate_infor=LearningAutomata(pop(i1).Position,idx_cat,newcost,lastcost,cate_infor);
                        end
                        
                        % Update Personal Best
                        if Dominates(pop(i1),pop(i1).Best)
                            pop(i1).Best.Position=pop(i1).Position;
                            pop(i1).Best.Cost=pop(i1).Cost;
                            pop(i1).Best.Out=pop(i1).Out; 
                            
                        elseif Dominates(pop(i1).Best,pop(i1))
                            % Do Nothing
                        else
                            if rand<0.5
                                pop(i1).Best.Position=pop(i1).Position;
                                pop(i1).Best.Cost=pop(i1).Cost;
                                pop(i1).Best.Out=pop(i1).Out; 
                            end
                        end                
                        m=m+1;
                    else
                        m=Ns;
                    end
                end % end swim
                
            end  % end chemotaxis of each bacterium
        
            % updating the second subpopulation-star topology----------------------
            hubind_temp=round(rand(1,1)*SubPoP(1));
            hubind=hubind_temp+SubPoP(1);
        
            for i2=SubPoP(1)+1:nPoP
                leader= rep([rep.ReporPopInd]==randperm(min(numel(rep),3),1)); 
                lastcost=pop(i2).Cost;
                
                pop(i2).Delta=2*rand(1,nVar)-1;         % random walk
                pop(i2).Tumble=pop(i2).Delta/sqrt(pop(i2).Delta*(pop(i2).Delta)');
                pop(i2).StepSize=Cmin+exp(-cc*(it/MaxIt)^2)*(Cmax-Cmin);
                
                tempmbest2=struct;
                tempmbest2(1).Position=pop(i2).Best.Position;
                tempmbest2(2).Position=pop(hubind).Best.Position;
                
                tempmbest2(1).Cost=pop(i2).Best.Cost;
                tempmbest2(2).Cost=pop(hubind).Best.Cost;
                
                tempmbest2=DetermineDomination(tempmbest2);
                tempzero=find([tempmbest2.IsDominated]==0);
                hubbest=tempmbest2(tempzero(1));
                
                % update numeric part
                if ~isempty(idx_num)
                    pop(i2).Position(:,idx_num)=pop(i2).Position(:,idx_num)+pop(i2).StepSize*(...
                        c1*rand.*(hubbest.Position(:,idx_num)-pop(i2).Position(:,idx_num))+...
                        c2*rand.*(leader.Position(:,idx_num)-pop(i2).Position(:,idx_num))+pop(i2).Tumble(:,idx_num));
                end
        
                % update categorical part
                if ~isempty(idx_cat)
                    pop(i1).Position(randperm(cn,1),idx_cat(j))=cate_infor(j).avaliable_total(betterProb(1));
                end
                
                % position limitations---------------------------------------------
                pop(i2) = boundarylimitations(pop(i2),idx_num,idx_cat,VarMax,VarMin);
                
                % compute fitness--------------------------------------------------
                [pop(i2).Cost, pop(i2).Out]=CostFunction(pop(i2).Position);    % calculate cost        
                newcost=pop(i2).Cost;
                
                % Learning Automata strategy to categorical attributes
                if ~isempty(idx_cat)
                    cate_infor=LearningAutomata(pop(i2).Position,idx_cat,newcost,lastcost,cate_infor);
                end
                
                % Update Personal Best
                if Dominates(pop(i2),pop(i2).Best)
                    pop(i2).Best.Position=pop(i2).Position;
                    pop(i2).Best.Cost=pop(i2).Cost;
                    pop(i2).Best.Out=pop(i2).Out; 
                    
                elseif Dominates(pop(i2).Best,pop(i2))
                    % Do Nothing
                else
                    if rand<0.5
                        pop(i2).Best.Position=pop(i2).Position;
                        pop(i2).Best.Cost=pop(i2).Cost;
                        pop(i2).Best.Out=pop(i2).Out; 
                    end
                end
                        
                % swim------------------------------------------------------
                m=0;
                while m<Ns
                    
                    if newcost <= lastcost             % compare cost to find better cost                
                       lastcost=newcost;
                       
                       tempmbest2(1).Position=pop(i2).Best.Position;
                       tempmbest2(2).Position=pop(hubind).Best.Position;
        
                       tempmbest2(1).Cost=pop(i2).Best.Cost;
                       tempmbest2(2).Cost=pop(hubind).Best.Cost;
        
                       tempmbest2=DetermineDomination(tempmbest2);
                       tempzero=find([tempmbest2.IsDominated]==0);
                       hubbest=tempmbest2(tempzero(1));
                       
                        % update numeric part
                        if ~isempty(idx_num)
                            pop(i2).Position(:,idx_num)=pop(i2).Position(:,idx_num)+pop(i2).StepSize*(...
                                c1*rand.*(hubbest.Position(:,idx_num)-pop(i2).Position(:,idx_num))+...
                                c2*rand.*(leader.Position(:,idx_num)-pop(i2).Position(:,idx_num)));
                        end
        
                        % update categorical part
                        if ~isempty(idx_cat)
                            pop(i1).Position(randperm(cn,1),idx_cat(j))=cate_infor(j).avaliable_total(betterProb(1));
                        end
        
                        % position limitations-------------------------------------
                        pop(i2) = boundarylimitations(pop(i2),idx_num,idx_cat,VarMax,VarMin);
                        
                        % compute fitness------------------------------------------
                        [pop(i2).Cost, pop(i2).Out]=CostFunction(pop(i2).Position);    % calculate cost                
                        newcost=pop(i2).Cost;
                        
                        % Learning Automata strategy to categorical attributes
                        if ~isempty(idx_cat)
                            cate_infor=LearningAutomata(pop(i2).Position,idx_cat,newcost,lastcost,cate_infor);
                        end
                        
                        % Update Personal Best
                        if Dominates(pop(i2),pop(i2).Best)
                            pop(i2).Best.Position=pop(i2).Position;
                            pop(i2).Best.Cost=pop(i2).Cost;
                            pop(i2).Best.Out=pop(i2).Out;
                            
                        elseif Dominates(pop(i2).Best,pop(i2))
                            % Do Nothing
                        else
                            if rand<0.5
                                pop(i2).Best.Position=pop(i2).Position;
                                pop(i2).Best.Cost=pop(i2).Cost;
                                pop(i2).Best.Out=pop(i2).Out; 
                            end
                        end
                                      
                        m=m+1;
                    else
                        m=Ns;
                    end
                end % end swim
            end % end the second subpopulation
            
            % get and update rep --------------------------------------------------
            rep=getRep(pop,rep,nRep,nGrid,alpha);
            
            RepCDInd=RepCrowdingDistance(rep);  % crowding distance
            
            [RepSOInd,total_SO]=getSO(rep);  % origin-nondominated solution distance
            
            if total_SOH==total_SO
                count_rep=count_rep+1;
            else
                count_rep=0;
                total_SOH=total_SO;
            end
            
            TempRepInd=RepCDInd'+RepSOInd;
            
            [~,RepInd]=sort(TempRepInd); % ascend
            
            for s=1:numel(rep)
                rep(s).ReporPopInd=RepInd(s);
            end

            % reproduction-elimination --------------------------------------------------------
            if mod(it,Fre)==0                % perform reproduction every Fre times
                pop=Reproduction_elimination(rep,pop,nPoP,Sr);
            end
            
            % elimination-dispersal -----------------------------------------------
            if mod(it,Fmig)==0 
                for i=1:nPoP
                    if rand<Pmig
                        for j=1:size(idx_num,2)
                            pop(i).Position(:,idx_num(j))=VarMin(1,idx_num(j))+(VarMax(1,idx_num(j))-VarMin(1,idx_num(j))).*rand(cn,1);
                        end
                        for j=1:size(idx_cat,2)        
                            for jj=1:cn
                                tempindex=randperm(size(data,1),1);
                                pop(i).Position(jj,idx_cat(j))=data(tempindex,idx_cat(j));
                            end
                        end    
                    end
                end
            end
            
            % deal with stagnation-------------------------------------------------    
            if count_rep>=1
                rep = Perturbation_rep(rep,CostFunction,VarMax,VarMin,idx_num,idx_cat,cn,data);
            end  
            
            % get and update rep --------------------------------------------------
            rep=getRep(pop,rep,nRep,nGrid,alpha);
                
            % Show Iteration Information
            disp(['TestProblem ' num2str(TP) ',Run ' num2str(run) '-Iteration '  num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
                    
        end % end Iteration
        
        %% obtain the results
        FianlOut=FianlSolution(rep);
        
        temp_Acc=FianlOut.Acc;
        temp_ARI=FianlOut.ARI;
        temp_Pre=FianlOut.Pre;
        temp_Fm=FianlOut.Fm;

        rep_all_runs{run, 1}=[rep.Cost];

        Acc_all_runs{run, 1} = temp_Acc;
        ARI_all_runs{run, 1} = temp_ARI;
        Pre_all_runs{run, 1} = temp_Pre;
        Fm_all_runs{run, 1} = temp_Fm;
                      
    end % end TotalRun/ parfor 

    %% obtain the evalution index results
    Acc_index=zeros(1,4);
    Acc_index(1)=mean(cell2mat(Acc_all_runs)); % average value
    [Acc_index(2),accindex]=max(cell2mat(Acc_all_runs)); % maximal value
    Acc_index(3)=min(cell2mat(Acc_all_runs)); % minimal value
    Acc_index(4)=std(cell2mat(Acc_all_runs)); % standard variance
    
    ARI_Index=zeros(1,4);
    ARI_Index(1)=mean(cell2mat(ARI_all_runs));
    ARI_Index(2)=max(cell2mat(ARI_all_runs));
    ARI_Index(3)=min(cell2mat(ARI_all_runs));
    ARI_Index(4)=std(cell2mat(ARI_all_runs));
    
    Pre_index=zeros(1,4);
    Pre_index(1)=mean(cell2mat(Pre_all_runs));
    Pre_index(2)=max(cell2mat(Pre_all_runs));
    Pre_index(3)=min(cell2mat(Pre_all_runs));
    Pre_index(4)=std(cell2mat(Pre_all_runs));
    
    Fm_Index=zeros(1,4);
    Fm_Index(1)=mean(cell2mat(Fm_all_runs));
    Fm_Index(2)=max(cell2mat(Fm_all_runs));
    Fm_Index(3)=min(cell2mat(Fm_all_runs));
    Fm_Index(4)=std(cell2mat(Fm_all_runs));
    
    %% save Resluts of each runtime
    columns1={'Accuracy','ARI','Precision','Fmeasure'};
    Runtimeresults=table(cell2mat(Acc_all_runs),cell2mat(ARI_all_runs),cell2mat(Pre_all_runs), cell2mat(Fm_all_runs),'VariableNames', columns1);
    writetable(Runtimeresults, [dataname  '_' datestr(now,'yyyy.mm.dd HH.MM') '_' 'Runtime_results.csv']);

    %% save final Resluts
    index={'mean';'max';'min';'std'};
    columns2={'index','Accuracy','ARI','Precision','Fmeasure'};
    finalresults=table(index, Acc_index',ARI_Index',Pre_index',Fm_Index','VariableNames',columns2);
    writetable(finalresults, [dataname '_' datestr(now,'yyyy.mm.dd HH.MM') '_' 'final_results.csv']);
    %% save all reps
    maxLength = max(cellfun(@length, rep_all_runs)); % Find the longest element length
    
    % Fill the short line to match the longest line length
    for i = 1:length(rep_all_runs)
        currentLength = length(rep_all_runs{i});
        if currentLength < maxLength
            rep_all_runs{i} = [rep_all_runs{i}, [1;1].*NaN(1, maxLength - currentLength)];
        end
    end
    
    rep_matrix = cell2mat(rep_all_runs);
    
    writematrix(rep_matrix, [dataname '_rep_output.csv']);
    
    %% save the best rep_Cost which accindex==1        
    final_rep_Cost=rep_all_runs{accindex};
    save(['BFO_' dataname '_' 'final_rep_Cost'],'final_rep_Cost');

end % end the test problem
