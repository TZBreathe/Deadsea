%% ARIMA model for capacity prediction

% Load data
clear;
Br{1}=gdFun.Load_Breathe_Run(507);

for i=1:length(Br)
    
    % create list of capacities to map to
    CapList{i}(1)=Br{i}.RunData.cycleTable{1,'ahDchrge'};
    for j=2:height(Br{i}.RunData.cycleTable)      
        CapList{i}(j)=Br{i}.RunData.cycleTable{j,'ahDchrge'};
        if CapList{i}(j)<CapList{i}(j-1)/1.5    % replace oddities in data where capacity suddenly drops for certain cycles
            CapList{i}(j)= CapList{i}(j-1);
            outlierCap{i}(j)=j;
        end
    end
end

CapList{1}(1)=CapList{1}(2); % temp fix bad initial data point
y=CapList{1}(:);

%ARIMA model
Mdl=arima(5,2,10);
Idx_est=100;
EstMdl=estimate(Mdl,y(1:Idx_est));

yf=forecast(EstMdl,196-Idx_est,y(1:Idx_est));
cycl=1:length(y);
figure;
plot(cycl,y,"b","LineWidth",2);
hold on
plot(cycl(Idx_est+1:end),yf,"r--","LineWidth",2);