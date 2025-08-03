
function [dataname,data,output,inputType,idx_num,new_idx_num,idx_cat,cn] = TestingDataSets(TP)

%% continuous data sets-----------------------------------------------------

% breast cancer coimbra
if TP==1
x=csvread('testing data\breast cancer coimbra.csv',1,0);

data = x(2:end,1:end-1);      
output = x(2:end,end);        

inputType = x(1,1:end-1);              % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp';

dataname='breast cancer coimbra';

cn=2;
end

% Ecoli
if TP==2
x=csvread('testing data\Ecoli.csv',1,0);

data = x(2:end,1:end-1);      %(last column, the output, is left out of the clustering)
output = x(2:end,end);        % the output (2 classes)

inputType = x(1,1:end-1);     % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 
 
dataname='Ecoli';

cn=8;
end

% glass
if TP==3
x=csvread('testing data\glass.csv',1,0);

data = x(2:end,1:end-1);      
output = x(2:end,end);        

inputType = x(1,1:end-1);              % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

dataname='glass';

cn=6;
end

% liver disorders
if TP==4
x=csvread('testing data\liver disorders.csv',1,0);

data = x(2:end,1:end-1);      
output = x(2:end,end);        

inputType = x(1,1:end-1);             % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 

dataname='Liver Disorders';

cn=2;
end

%% categorical data set
% vote
if TP==5
x=csvread('testing data\vote.csv',1,0);

data = x(2:end,1:end-1);      %(last column, the output, is left out of the clustering)
output = x(2:end,end);        % the output (2 classes)

inputType = x(1,1:end-1);     % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 
 
dataname='vote';

cn=2;
end

% zoo
if TP==6
x=csvread('testing data\zoo.csv',1,0);

data = x(2:end,1:end-1);      %(last column, the output, is left out of the clustering)
output = x(2:end,end);        % the output (2 classes)

inputType = x(1,1:end-1);     % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 
 
dataname='zoo';

cn=7;
end

%% mixed data sets----------------------------------------------------------
% australian
if TP==7
x=csvread('testing data\australian.csv',1,0);

data = x(2:end,1:end-1);      %(last column, the output, is left out of the clustering)
output = x(2:end,end);        % the output (2 classes)

inputType = x(1,1:end-1);     % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 

dataname='australian';

cn=2;
end

% german
if TP==8
x=csvread('testing data\german.csv',1,0);

data = x(2:end,1:end-1);      %(last column, the output, is left out of the clustering)
output = x(2:end,end);        % the output (2 classes)

inputType = x(1,1:end-1);     % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 

dataname='german';

cn=2;
end

% heart failure clinical records
if TP==9
x=csvread('testing data\heart failure clinical records.csv',1,0);

data = x(2:end,1:end-1);      %(last column, the output, is left out of the clustering)
output = x(2:end,end);        % the output (2 classes)

inputType = x(1,1:end-1);     % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 

dataname='heart failure clinical records';

cn=2;
end

% heart
if TP==10
x=csvread('testing data\heart.csv',1,0);

data = x(2:end,1:end-1);      %(last column, the output, is left out of the clustering)
output = x(2:end,end);        % the output (2 classes)

inputType = x(1,1:end-1);     % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 

dataname='Heart';

cn=2;
end

% japan_heart_attack_dataset -preprocessing  // kaggle
if TP==11
x=csvread('testing data\Heart Attack in Japan.csv',1,0);

data = x(2:end,1:end-1);      %(last column, the output, is left out of the clustering)
output = x(2:end,end);        % the output (2 classes)

inputType = x(1,1:end-1);     % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 

dataname='Heart Attack in Japan';

cn=2;
end

% heart_attack_youth_vs_adult_pre
if TP==12
x=csvread('testing data\Heart Attack in America(State).csv',1,0);

data = x(2:end,1:end-1);      %(last column, the output, is left out of the clustering)
output = x(2:end,end);        % the output (2 classes)

inputType = x(1,1:end-1);     % 1 for categorical, 0 for numeric
idx_cat = find(inputType);
idx_num = find(~inputType);
new_idx_num=[idx_num size(data,2)+1];

temp=mapminmax(data(:,idx_num)',0,1); % stanadralization
data(:,idx_num)=temp'; 

dataname='Heart Attack in America(State)';

cn=2;
end

end
