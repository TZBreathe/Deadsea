%% starter code for extracting features for degradation regression model

% Varta Client data
addpath(genpath('G:\09_Research_Development_Projects\00_Degradation data sets'));
filename=cell(6,1);
filename={'Varta188_07C_415V','Varta189_07C_415V','Varta190_07C_415V','Varta191_05C_415V','Varta192_05C_415V','Varta193_05C_415V'};
cyclelife=[316 369 435 1672 1254 1657];  % Defined as cycle to 85% BoL
MaxCap=3;
% variance of dQ(V)

for i=1:length(filename)
    
    Data=struct2cell(load(filename{i}));
    ind_10=find([Data{1, 1}.RunData.dataTable{:,'cycleNumber'}==5]&(Data{1, 1}.RunData.dataTable{:,'currCell'}<0));
    ind_50=find([Data{1, 1}.RunData.dataTable{:,'cycleNumber'}==51]&(Data{1, 1}.RunData.dataTable{:,'currCell'}<0));
     ind_100=find([Data{1, 1}.RunData.dataTable{:,'cycleNumber'}==450]&(Data{1, 1}.RunData.dataTable{:,'currCell'}<0));
     
    VQ_10=[Data{1, 1}.RunData.dataTable{ind_10,'voltCell'} Data{1, 1}.RunData.dataTable{ind_10,'ahDchge'}-Data{1, 1}.RunData.dataTable{ind_10(1),'ahDchge'}];
    VQ_10(:,2)=VQ_10(:,2)./MaxCap; %Normalise capacity, make sense or not?
    V_interp=linspace(VQ_10(1,1),VQ_10(end,1),100);
    Q_interp=interp1(VQ_10(:,1),VQ_10(:,2),V_interp);
    VQ_10=[V_interp; Q_interp];
    
    VQ_50=[Data{1, 1}.RunData.dataTable{ind_50,'voltCell'} Data{1, 1}.RunData.dataTable{ind_50,'ahDchge'}-Data{1, 1}.RunData.dataTable{ind_50(1),'ahDchge'}];
    VQ_50(:,2)=VQ_50(:,2)./MaxCap;
    V_interp=linspace(VQ_50(1,1),VQ_50(end,1),100);
    Q_interp=interp1(VQ_50(:,1),VQ_50(:,2),V_interp);
    VQ_50=[V_interp; Q_interp];
    
     dQ=abs(VQ_10(2,:)-VQ_50(2,:));
     VarQ(i)=var(dQ);
%      AvgQ(i)=mean(dQ);
%      MaxQ(i)=max(dQ);
     
     EndCap(i)=(Data{1, 1}.RunData.dataTable{ind_100(end),'ahDchge'}-Data{1, 1}.RunData.dataTable{ind_100(1),'ahDchge'})./MaxCap; 
     
     %     plot(VQ_10(1,:),VQ_10(2,:),VQ_50(1,:),VQ_50(2,:));
%        plot(dQ,VQ_4(1,:));

end
     

%  plot(EndCap, VarQ,'O');
% figure
% plot(cyclelife, AvgQ,'O');
% figure
% plot(cyclelife, MaxQ,'O');

%%

