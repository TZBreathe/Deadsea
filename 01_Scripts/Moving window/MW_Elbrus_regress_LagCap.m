%% Moving window prediction preliminary regression 

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

CapList={};
dVarQ={};
Cap_future={};
outlierCap={};
dNegPot={};
dChTime={};
% outlierNegPot={};
% outlierdQ={};
% outlierdChTime={};


%% Extract features - List of capacities & Var(dQ)


for i=1:length(Br)
    
    % create list of capacities to map to
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
Future=50;
Past=40;
%%
%var(dQ)
for i=1:length(Br)

    ii=1; %index for storing data points each cycle
    
    for j=Now:(length(CapList{i})-Now)  % index for 'now', begin at cycle 50 
    
     Cap_future{i}(ii)=CapList{i}(j+Future)/MaxCap; % Capacities of 'future' (j+50) cycles normalised against MaxCap
     Cap_lag0{i}(ii)=CapList{i}(j)/MaxCap; % Lagged capacities used as training input
     Cap_lag5{i}(ii)=CapList{i}(j-5)/MaxCap;
     Cap_lag10{i}(ii)=CapList{i}(j-10)/MaxCap;
%    dCap{i}(ii)=CapList{i}(j)-CapList{i}(j+50); % Difference in capacity between now and 50 cycles later
            
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
    [~,ia,ic]=unique(VQ_now_tmp(:,1));   %remove repeated V values
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

clear dQ VQ_early_tmp  VQ_early V_interp Q_interp VQ_now_tmp VQ_now ia ic
   %% Avg negPot. Can use charge and/or discharge
   
 
for i=1:length(Br)
    ii=1;  
    for j=Now:(length(CapList{i})-Now)    
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
    for j=Now:1:(length(CapList{i})-Now)    

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

%% lasso regression 
clear Br
idx=[1 2 3 4 5 6 7];
xdVarQ=[];xdNegPot=[];xdChTime=[];xCap_lag0=[];xCap_lag5=[];xCap_lag10=[];yCap=[];
%Arrange training dat
for i =1:length(idx)
    xdVarQ=cat(1,xdVarQ,dVarQ{idx(i)}');
    xdNegPot=cat(1,xdNegPot,dNegPot{idx(i)}');
    xdChTime=cat(1,xdChTime,dChTime{idx(i)}');
    xCap_lag0=cat(1,xCap_lag0,Cap_lag0{idx(i)}');
    xCap_lag5=cat(1,xCap_lag5,Cap_lag5{idx(i)}');
    xCap_lag10=cat(1,xCap_lag10,Cap_lag10{idx(i)}');
    yCap=cat(1,yCap,Cap_future{idx(i)}');
end
[xdVarQ,cdVarQ,sdVarQ]=normalize(xdVarQ);
[xdNegPot,cdNegPot,sdNegPot]=normalize(xdNegPot);
[xdChTime,cdChTime,sdChTime]=normalize(xdChTime);
[xCap_lag0,cCap_lag0, sCap_lag0]=normalize(xCap_lag0);
[xCap_lag5,cCap_lag5, sCap_lag5]=normalize(xCap_lag5);
[xCap_lag10,cCap_lag10, sCap_lag10]=normalize(xCap_lag10);
X=[xdVarQ xdNegPot xdChTime xCap_lag0 xCap_lag5 xCap_lag10];

% fit lasso regression
[b1,FitInfo]=lasso(X,yCap,'Alpha',0.7,'Intercept',true,'CV',10);
idxLambda1SE = FitInfo.Index1SE;
coef = b1(:,idxLambda1SE);
coef0 = FitInfo.Intercept(idxLambda1SE);

yhat = X*coef + coef0;
hold on;

scatter(yCap,yhat);plot(yCap,yCap);
xlabel('True capacity');ylabel('Fitted capacity');
hold off;
mean(abs(rmmissing(yhat-yCap)))


% scatter(yCap,yhat-yCap);

%% Test regression

% define test data set
j=8;
xdVarQ_v=normalize([dVarQ{j}]','center',cdVarQ,'scale',sdVarQ);
xdNegPot_v=normalize(dNegPot{j}','center',cdNegPot,'scale',sdNegPot);
xdChTime_v=normalize(dChTime{j}','center',cdChTime,'scale',sdChTime);
xCap_lag0_v=normalize(Cap_lag0{j}','center',cCap_lag0,'scale',sCap_lag0);
xCap_lag5_v=normalize(Cap_lag5{j}','center',cCap_lag5,'scale',sCap_lag5);
xCap_lag10_v=normalize(Cap_lag10{j}','center',cCap_lag10,'scale',sCap_lag10);
yCap_v=[Cap_future{j}'];
X_v=[xdVarQ_v xdNegPot_v xdChTime_v xCap_lag0_v xCap_lag5_v xCap_lag10_v];


yhat_v = X_v*coef + coef0;

hold on;
scatter(yCap_v,yhat_v);plot(yCap_v,yCap_v);
xlabel('True capacity');ylabel('Fitted capacity');
% xlim([min(yCap_v) max(yCap_v)]); ylim([min(yCap_v) max(yCap_v)]);
hold off;
mean(abs(rmmissing(yhat_v-yCap_v)))

% scatter(yCap_v,yhat_v-yCap_v);

% 

