%% ARIMAX model for capacity prediction

% Load data
clear;
Br{1}=gdFun.Load_Breathe_Run(506);
Br{2}=gdFun.Load_Breathe_Run(507);
Br{3}=gdFun.Load_Breathe_Run(508);
Br{4}=gdFun.Load_Breathe_Run(509);

runNumbers=[767 766 846]; 
Br{5}=appendDegData(runNumbers);
runNumbers=[765 764 845]; 
Br{6}=appendDegData(runNumbers);
runNumbers=[771 769 848]; 
Br{7}=appendDegData(runNumbers);
runNumbers=[776 773 850]; 
Br{8}=appendDegData(runNumbers);

% runNumbers=[763 844]; 
% Br{9}=appendDegData(runNumbers);
runNumbers=[1101 1104]; 
Br{9}=appendDegData(runNumbers);
runNumbers=[1107 1106]; 
Br{10}=appendDegData(runNumbers);
Br{11}=gdFun.Load_Breathe_Run(1105);
Br{12}=gdFun.Load_Breathe_Run(1102);
Br{13}=gdFun.Load_Breathe_Run(1103);
Br{14}=gdFun.Load_Breathe_Run(1108);

MaxCap=4.9;
n=3; %skip every n data points

%%
CapList={};
dVarQ={};
Cap_future={};
outlierCap={};
dChTime={};
n=3; %skip every n data points
y=[];
yS=[];
for i=1:length(Br)
    
    % create list of capacities to map to
    CapList{i}(1)=Br{i}.RunData.cycleTable{1,'ahDchrge'};
    for j=2:height(Br{i}.RunData.cycleTable)/n      
        CapList{i}(j)=Br{i}.RunData.cycleTable{n*j,'ahDchrge'};
        if CapList{i}(j)<CapList{i}(j-1)/1.5    % replace oddities in data where capacity suddenly drops for certain cycles
            CapList{i}(j)= CapList{i}(j-1);
            outlierCap{i}(j)=j;
        end
    end
   
    CapList{i}(1)=CapList{i}(2); % temp fix bad initial data point 
    y{i}=CapList{i}(:)/MaxCap;   
    yS{i}=smooth(y{i},4); %smoothing
end

%% ARIMA model

Idx_est=36;
Mdl=arima(3,2,2);


for i=1:10
[yP{i}]=polyfit(1:Idx_est,y{i}(1:Idx_est),1);
EstMdl(i)=estimate(Mdl,yS{i}(1:Idx_est));
end

for i=1:10
yf{i}=forecast(EstMdl(i),length(y{i})-Idx_est,yS{i}(1:Idx_est));
yP{i}=polyval(yP{i},1:length(y{i}));
cycl{i}=1:length(y{i});
end

subplot(4,2,1);
hold on;
plot(cycl{1}*n,y{1},"b","LineWidth",2);
plot(cycl{1}(Idx_est+1:end)*n,yf{1},"r--","LineWidth",2);
plot(cycl{1}*n,yP{1},"g--","LineWidth",2);
xlabel('Cycles');ylabel('%Capacity');
legend('Actual','ARIMA','Linear extrapolation','location','best');

for i=1:5
subplot(3,2,i+1);
hold on;
plot(cycl{2*i}*n,y{2*i},"b","LineWidth",2);
plot(cycl{2*i}(Idx_est+1:end)*n,yf{2*i},"r--","LineWidth",2);
plot(cycl{2*i}*n,yP{2*i},"g--","LineWidth",2);
xlabel('Cycles');ylabel('%Capacity');

end

hold off;

