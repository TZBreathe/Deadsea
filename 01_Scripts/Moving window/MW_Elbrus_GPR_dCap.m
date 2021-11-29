%% Moving window GPR Elbrus data
%Load data
clear;
MaxCap=4.9;
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
% Br{16}=gdFun.Load_Breathe_Run(1109);


CapList={};
dVarQ={};
Cap_future={};
outlierCap={};
dNegPot={};
dChTime={};
% outlierNegPot={};
% outlierdQ={};
% outlierdChTime={};

%% List of capacities & Var(dQ) and initialise


for i=1:length(Br)
    
    % create list of capacities 
    CapList{i}(1)=Br{i}.RunData.cycleTable{1,'ahDchrge'};
    for j=2:height(Br{i}.RunData.cycleTable)      
        CapList{i}(j)=Br{i}.RunData.cycleTable{j,'ahDchrge'};
        if CapList{i}(j)<CapList{i}(j-1)/1.5    % replace oddities in data where capacity suddenly drops for certain cycles
            CapList{i}(j)= CapList{i}(j-1);
%             outlierCap{i}(j)=j;
        end
    end
end

Now=50;
Future=40;
Past=40;

%%
%var(dQ)
for i=1:length(Br)

    ii=1; %index for storing data points each cycle
    
    for j=Now:(length(CapList{i})-Future)  % index for 'now', begin at cycle 50 
    
%      Cap_future{i}(ii)=CapList{i}(j+Future)/MaxCap; % Capacities of 'future' (j+50) cycles normalised against MaxCap
      Cap_lag0{i}(ii)=CapList{i}(j)/MaxCap; % Lagged capacities used as training input
%      Cap_lag5{i}(ii)=CapList{i}(j-5)/MaxCap;
%      Cap_lag10{i}(ii)=CapList{i}(j-10)/MaxCap;
      dCap{i}(ii)=(CapList{i}(j)-CapList{i}(j+Future))/MaxCap; % Difference in capacity between now and 50 cycles later
      dCap_lag10{i}(ii)=abs(CapList{i}(j)-CapList{i}(j-10))/MaxCap;    
      dCap_lag40{i}(ii)=abs(CapList{i}(j)-CapList{i}(j-40))/MaxCap; 
      
      
    %Find early (j-40) and now (j) indexes for voltage and capacity
    idx_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j-Past]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
    idx_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&(Br{i}.RunData.dataTable{:,'currCell'}<0)); 
    
    %some datasets have missing discharge sequences, drop these
    if length(idx_early)<10 | length(idx_now)<10
      dVarQ{i}(ii)=NaN;
      ii=ii+1;
      continue;
    end 
    % Get capacity and voltage for 'early' (j-40) cycle, then reshape to
    % same length
    VQ_early_tmp=[Br{i}.RunData.dataTable{idx_early,'voltCell'} Br{i}.RunData.dataTable{idx_early,'ahDchrge'}];
    [~,ia,~]=unique(VQ_early_tmp(:,1));   %remove repeated V values so interp1 could work
    VQ_early=VQ_early_tmp(ia,:);
    VQ_early(:,2)=VQ_early(:,2)./MaxCap; %Normalise capacity with assumed Max value for all cells of same type
    %Interpolate into same length for all 
    V_interp=linspace(VQ_early(1,1),VQ_early(end,1),1000);
    Q_interp=interp1(VQ_early(:,1),VQ_early(:,2),V_interp);
    VQ_early=[V_interp; Q_interp];
    % Get capacity and voltage for 'now'
    VQ_now_tmp=[Br{i}.RunData.dataTable{idx_now,'voltCell'} Br{i}.RunData.dataTable{idx_now,'ahDchrge'}];
    [~,ia,~]=unique(VQ_now_tmp(:,1));   %remove repeated V values
    VQ_now=VQ_now_tmp(ia,:);
    VQ_now(:,2)=VQ_now(:,2)./MaxCap;
    V_interp=linspace(VQ_now(1,1),VQ_now(end,1),1000);
    Q_interp=interp1(VQ_now(:,1),VQ_now(:,2),V_interp);
    VQ_now=[V_interp; Q_interp];
    
   %compute var(dQ)
    dQ=abs(VQ_early(2,:)-VQ_now(2,:));
    dVarQ{i}(ii)=(var(dQ));
    
    ii=ii+1;
    
    end
    
    %Some cycle data contains oddities resulting in outliers in var(dQ),
    %find these outliers & their indexes
    %replace outlier with NaN so they don't affect regression. 
    [dVarQ{i},outlierdQ_tmp]=filloutliers(dVarQ{i},NaN); 
%     outlierdQ{i}=outlierdQ_tmp;
     
end

clear dQ VQ_early_tmp  VQ_early V_interp Q_interp VQ_now_tmp VQ_now ia 

   %% Avg negPot. Can use charge and/or discharge
 
for i=1:length(Br)
    ii=1;  
    for j=Now:(length(CapList{i})-Future)    
    k=1; % index 1-3 for storing avearge over 3 cycles
        for l=(j-Past-2):(j-Past) %Average over 3 cycles for stable readings
        idx_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==l]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
            if isempty(idx_early) 
             continue;
            end
         NegPot(k)=mean(Br{i}.RunData.dataTable{idx_early,'surpotNeg'});
         k=k+1;
        end
     NegPotAvg_early=mean(NegPot);
     
    k=1; 
        for l=(j-2):j  %Average over 3 cycles for stable readings
         idx_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==l]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
            if isempty(idx_now) 
            continue;
            end
          NegPot(k)=mean(Br{i}.RunData.dataTable{idx_now,'surpotNeg'});
            k=k+1;
        end
    NegPotAvg_now=mean(NegPot);
    
    dNegPot{i}(ii)=NegPotAvg_early-NegPotAvg_now;
    ii=ii+1;
    end
    [dNegPot{i},outlierNegPot_tmp]=filloutliers(dNegPot{i},NaN); %set outliers to NaN to not affect regression
%     outlierNegPot{i}=outlierNegPot_tmp;
end

clear NegPotAvg_early NegPotAvg_now NegPot

%%
% Totalcharge time
for i=1:length(Br)
 
    ii=1;
    for j=Now:1:(length(CapList{i})-Future)    

    idx_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j-Past]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
    idx_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
    
      if isempty(idx_early) | isempty(idx_now)
       dChTime{i}(ii)=0;
       ii=ii+1;
      continue;
      end 
    
    ChTime_early=Br{i}.RunData.dataTable{idx_early(end),'timeTest'}-Br{i}.RunData.dataTable{idx_early(1),'timeTest'};
    Chtime_now=Br{i}.RunData.dataTable{idx_now(end),'timeTest'}-Br{i}.RunData.dataTable{idx_now(1),'timeTest'};
    
    dChTime{i}(ii)=Chtime_now/ChTime_early;
    ii=ii+1;    
    end

    [dChTime{i},outlierdChTime_tmp]=filloutliers(dChTime{i},NaN);
%     outlierdChTime{i}=outlierdChTime_tmp;  
end
clear  ChTime_early Chtime_now

%% check outliers

% Cap_outlier=find(outlierCap{1});
% dQ_outlier=find(outlierdQ{1});
% ChTime_outlier=find(outlierdChTime{1});
% NegPot_outlier=find(outlierNegPot{1});

%% GPR 

% clear Br
idx=[1 3 4 6 7 8 9 10 14 12 13];
xdVarQ=[];xdNegPot=[];xdChTime=[];xCap_lag0=[];xdCap_lag10=[];xdCap_lag40=[];ydCap=[];

%Arrange training data
for i =1:length(idx)
    xdVarQ=cat(1,xdVarQ,dVarQ{idx(i)}');
    xdNegPot=cat(1,xdNegPot,dNegPot{idx(i)}');
    xdChTime=cat(1,xdChTime,dChTime{idx(i)}');
    xCap_lag0=cat(1,xCap_lag0,Cap_lag0{idx(i)}');
    xdCap_lag10=cat(1,xdCap_lag10,dCap_lag10{idx(i)}');
    xdCap_lag40=cat(1,xdCap_lag40,dCap_lag40{idx(i)}');
    ydCap=cat(1,ydCap,dCap{idx(i)}');
end

[X,C,S]=normalize([xdVarQ xdNegPot xdChTime xCap_lag0 xdCap_lag10 xdCap_lag40]);

%Fit GPR
gprMdl=fitrgp(X,ydCap,'Basis','linear','FitMethod','exact','PredictMethod','exact');

%Predict with training data
yhat = predict(gprMdl,X);

%% Test regression

xdVarQ_v=[];xdChTime_v=[];xCap_lag0_v=[];xdCap_lag10_v=[];xdCap_lag40_v=[]; ydCap_v=[];X_v=[];xdNegPot_v=[];
idx_v=[2 5 11];
for i =1:length(idx_v)
    xdVarQ_v=cat(1,xdVarQ_v,dVarQ{idx_v(i)}');
    xdNegPot_v=cat(1,xdNegPot_v,dNegPot{idx_v(i)}');
    xdChTime_v=cat(1,xdChTime_v,dChTime{idx_v(i)}');
    xCap_lag0_v=cat(1,xCap_lag0_v,Cap_lag0{idx_v(i)}');
    xdCap_lag10_v=cat(1,xdCap_lag10_v,dCap_lag10{idx_v(i)}');
    xdCap_lag40_v=cat(1,xdCap_lag40_v,dCap_lag40{idx_v(i)}');
    ydCap_v=cat(1, ydCap_v,dCap{idx_v(i)}');
end

X_v=normalize([xdVarQ_v xdNegPot_v xdChTime_v xCap_lag0_v xdCap_lag10_v xdCap_lag40_v],'center',C,'scale',S);

yhat_v = predict(gprMdl,X_v);


hold on;
scatter(ydCap*100,yhat*100);
xlabel('True capacity');ylabel('Fitted capacity');
% hold off;

scatter(ydCap_v*100,yhat_v*100);
xlabel('True capacity change %');ylabel('Predicted capacity change %');
plot(ydCap*100,ydCap*100,'red');
legend('Traning','Validation','True value','location','best');
hold off;
mean(abs(rmmissing(yhat-ydCap)))
mean(abs(rmmissing(yhat_v-ydCap_v)))



% 

