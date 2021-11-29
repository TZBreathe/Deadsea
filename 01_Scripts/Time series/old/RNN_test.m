%% ARIMAX model for capacity prediction

% Load data
clear;
% Br{1}=gdFun.Load_Breathe_Run(506);
 Br{1}=gdFun.Load_Breathe_Run(507);
% Br{3}=gdFun.Load_Breathe_Run(508);
% Br{4}=gdFun.Load_Breathe_Run(509);

MaxCap=4.9;
n=3; %skip every n data points

%%
CapList={};
dVarQ={};
Cap_future={};
outlierCap={};
dChTime={};

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

%% RNN model

Idx_est=40;
% XTrain=[y{1}(1:Idx_est-3)';y{1}(2:Idx_est-2)';y{1}(3:Idx_est-1)'];
XTrain=[yS{1}(1:Idx_est-2)';yS{1}(2:Idx_est-1)'];
YTrain=[yS{1}(3:Idx_est)'];

mu = mean(yS{1}(1:Idx_est));
sig = std(yS{1}(1:Idx_est));
XTrain=(XTrain-mu)/sig;
YTrain=(YTrain-mu)/sig;
XTest=[y{1}(Idx_est+1:end-2)';y{1}(Idx_est+2:end-1)'];
YTest=y{1}(Idx_est+3:end)';

numFeatures = 2;
numResponses = 1;
numHiddenUnits = 300;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',250, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.01, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress');

net = trainNetwork(XTrain,YTrain,layers,options);

YPred=predict(net, XTrain);
plot(YPred);hold on;plot(YTrain);
%% Test model
XTest=(XTest-mu)/sig;
net = predictAndUpdateState(net,XTrain);
[net,YPred] = predictAndUpdateState(net,[YTrain(end-1); YTrain(end)]);
[net,YPred(:,2)] = predictAndUpdateState(net,[YTrain(end);YPred(:,1)]);
numTimeStepsTest = numel(YTest);

for i = 3:numTimeStepsTest

    [net,YPred(:,i)] = predictAndUpdateState(net,[YPred(:,i-2);YPred(:,i-1)],'ExecutionEnvironment','cpu');

end

YPred = sig*YPred + mu;
%  
% 
plot(YPred);hold on;plot(YTest);


YPred2=predict(net, XTest);
YPred2=sig*YPred2 + mu;
plot(YPred2);hold on;plot(YTest);

%%






