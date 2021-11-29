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
   
    CapList{i}(1)=CapList{i}(3); % temp fix bad initial data point 
    y{i}=CapList{i}(:);
   yS{i}=smooth(y{i});
end


%% Exogenous features

% Var(dQ)
for i=1:length(Br)

    ii=2; %another index for storing data points each cycle
    
    for j=2:height(Br{i}.RunData.cycleTable)/n  % index for 'now', begin at cycle 5
    
%     Cap_future{i}(ii)=CapList{i}(j+50); % Capacities of 'future' (j+50) cycles
%       dCap{i}(ii)=CapList{i}(j)-CapList{i}(j+50); % Difference in capacity between now and 50 cycles later
         
    %Find early (j-5) and now (j) indexes for voltage and capacity
    ind_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==n*j-n]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
    ind_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==n*j]&(Br{i}.RunData.dataTable{:,'currCell'}<0)); 
    
    %some datasets have missing discharge sequences, drop these
    if length(ind_early)<10 | length(ind_now)<10
      VarQ{i}(ii)=NaN;
      ii=ii+1;
      continue;
    end 
    % Get capacity and voltage for 'early' (j-40) cycle, then reshape to
    % same length
    VQ_early_tmp=[Br{i}.RunData.dataTable{ind_early,'voltCell'} Br{i}.RunData.dataTable{ind_early,'ahDchrge'}];
    [~,ia,~]=unique(VQ_early_tmp(:,1));   %remove repeated V values so interp1 could work
    VQ_early=VQ_early_tmp(ia,:);
    VQ_early(:,2)=VQ_early(:,2)./MaxCap; %Normalise capacity with assumed Max value for all cells of same type
    %Interpolate into same length for all 
    V_interp=linspace(VQ_early(1,1),VQ_early(end,1),1000);
    Q_interp=interp1(VQ_early(:,1),VQ_early(:,2),V_interp);
    VQ_early=[V_interp; Q_interp];
    % Get capacity and voltage for 'now'
    VQ_now_tmp=[Br{i}.RunData.dataTable{ind_now,'voltCell'} Br{i}.RunData.dataTable{ind_now,'ahDchrge'}];
    [~,ia,ic]=unique(VQ_now_tmp(:,1));   %remove repeated V values
    VQ_now=VQ_now_tmp(ia,:);
    VQ_now(:,2)=VQ_now(:,2)./MaxCap;
    V_interp=linspace(VQ_now(1,1),VQ_now(end,1),1000);
    Q_interp=interp1(VQ_now(:,1),VQ_now(:,2),V_interp);
    VQ_now=[V_interp; Q_interp];
    
   %compute var(dQ)
    dQ=abs(VQ_early(2,:)-VQ_now(2,:));
    VarQ{i}(ii)=(var(dQ));
    
    ii=ii+1;
    
    end
     VarQ{i}(1)= VarQ{i}(2); %tmp fix
    %Some cycle data contains oddities resulting in outliers in var(dQ),
    %find these outliers & their indexes
    %replace outlier with NaN so they don't affect regression. 
    
    [VarQ{i},outlierdQ_tmp]=filloutliers(VarQ{i},'next'); 
    
     
end

%% Exogenous features

% NegPot
for i=1:length(Br)
    ii=2;  
    for j=2:height(Br{i}.RunData.cycleTable)/n
        ind_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==n*j-n]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
%             if isempty(ind_early) 
%              continue;
%             end
         NegPotAvg_early=mean(Br{i}.RunData.dataTable{ind_early,'surpotNeg'});
         ind_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==n*j]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
%             if isempty(ind_now) 
%             continue;
%             end
          NegPotAvg_now=mean(Br{i}.RunData.dataTable{ind_now,'surpotNeg'});
                     
    dNegPot{i}(ii)=NegPotAvg_early-NegPotAvg_now;
    ii=ii+1;
    end
    dNegPot{i}(1)=dNegPot{i}(1);
    [dNegPot{i},outlierNegPot_tmp]=filloutliers(dNegPot{i},'next'); %set outliers to NaN to not affect regression
    
end

%% Exogenous features

% ChTime
for i=1:length(Br)
 
    ii=2;
    for j=2:height(Br{i}.RunData.cycleTable)/n  

    ind_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==n*j-n]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
    ind_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==n*j]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
    
      if isempty(ind_early) | isempty(ind_now)
       dChTime{i}(ii)=0;
       ii=ii+1;
      continue;
      end 
    
    ChTime_early=Br{i}.RunData.dataTable{ind_early(end),'timeTest'}-Br{i}.RunData.dataTable{ind_early(1),'timeTest'};
    Chtime_now=Br{i}.RunData.dataTable{ind_now(end),'timeTest'}-Br{i}.RunData.dataTable{ind_now(1),'timeTest'};
    
    dChTime{i}(ii)=Chtime_now/ChTime_early;
    ii=ii+1;    
    end
    dChTime{i}(2)=dChTime{i}(1);
    [dChTime{i},outlierdChTime_tmp]=filloutliers(dChTime{i},'next');
     
end

%% examine and forcast exogenous features
% 

%normalise exogeneous variables
xVarQ={};
xdNegPot={};
xdChTime={};
X={};
X_f={};

Idx_est=30; %Set number of cycles for training. actual cycle no. = Idx_test*n

for i=1:3
xVarQ{i}=normalize([VarQ{i}']);
xdNegPot{i}=normalize([dNegPot{i}']);
xdChTime{i}=normalize([dChTime{i}']);


% forcast exgeneous variables
Mdl_VarQ=arima(2,2,2);
Mdl_exo=arima(3,1,3);

EstMdl_varQ=estimate(Mdl_exo,xVarQ{i});
f_varQ=forecast(EstMdl_varQ,length(xVarQ{i})-Idx_est,xVarQ{i}(1:Idx_est));

EstMdl_NegPot=estimate(Mdl_exo,xdNegPot{i});
f_NegPot=forecast(EstMdl_NegPot,length(xdNegPot{i})-Idx_est,xdNegPot{i}(1:Idx_est));

EstMdl_ChTime=estimate(Mdl_exo,xdChTime{i});
f_ChTime=forecast(EstMdl_ChTime,length(xdChTime{i})-Idx_est,xdChTime{i}(1:Idx_est));

% X{i}=[xVarQ{i} xdNegPot{i} xdChTime{i}];
% X_f{i}=[f_varQ f_NegPot f_ChTime];
X{i}=[xVarQ{i} xdNegPot{i}];
X_f{i}=[f_varQ f_NegPot];
end
% 
% cycl=1:length(xVarQ{i});
% figure;
% plot(cycl,xVarQ{i},"b","LineWidth",2);
% hold on;
% plot(cycl(Idx_est+1:end),f_varQ,"r--","LineWidth",2);
%  xlabel('Cycle');ylabel('Var(dQ)');
% hold off;

% figure;
% plot(cycl,xdNegPot{i},"b","LineWidth",2);
% hold on;
% plot(cycl(Idx_est+1:end),f_NegPot,"r--","LineWidth",2);
% xlabel('Cycle');ylabel('Difference in avearge NE potential');
% hold off;
% % % 
% figure;
% plot(cycl*3,xdChTime{i},"b","LineWidth",2);
% hold on;
% plot(cycl(Idx_est+1:end)*3,f_ChTime,"r--","LineWidth",2);
% xlabel('Cycle');ylabel('Charge time difference');
% hold off;
% % 

%% ARIMAX model


Mdl=arima(3,2,2);


EstMdl1=estimate(Mdl,yS{1}(1:Idx_est),'X',X{1});
EstMdl2=estimate(Mdl,yS{2}(1:Idx_est),'X',X{2});
EstMdl3=estimate(Mdl,yS{3}(1:Idx_est),'X',X{3});

yf1=forecast(EstMdl1,length(y{1})-Idx_est,yS{1}(1:Idx_est));
yf1_exo=forecast(EstMdl1,length(y{1})-Idx_est,yS{1}(1:Idx_est),'X0',X{1},'XF',X_f{1});

yf2=forecast(EstMdl2,length(y{2})-Idx_est,yS{2}(1:Idx_est));
yf2_exo=forecast(EstMdl2,length(y{2})-Idx_est,yS{2}(1:Idx_est),'X0',X{2},'XF',X_f{2});

yf3=forecast(EstMdl3,length(y{3})-Idx_est,yS{3}(1:Idx_est));
yf3_exo=forecast(EstMdl3,length(y{3})-Idx_est,yS{3}(1:Idx_est),'X0',X{3},'XF',X_f{3});
% yf4=forecast(EstMdl,length(y1)-Idx_est,[y1(1:Idx_est) y2(1:Idx_est) y3(1:Idx_est)]);

cycl=1:length(y{1});
figure;hold on;
plot(cycl*3,y{1},"b","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf1,"r--","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf1_exo,"g--","LineWidth",2);
xlabel('Cycles');ylabel('Capacity');
legend('Exp','ARIMA','ARIMAX');
hold off;

cycl=1:length(y{2});
figure;hold on;
plot(cycl*3,y{2},"b","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf2,"r--","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf2_exo,"g--","LineWidth",2);
xlabel('Cycles');ylabel('Capacity');
legend('Exp','ARIMA','ARIMAX');
hold off;

cycl=1:length(y{3});
figure;hold on;
plot(cycl*3,y{3},"b","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf3,"r--","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf3_exo,"g--","LineWidth",2);
xlabel('Cycles');ylabel('Capacity');
legend('Exp','ARIMA','ARIMAX');
hold off;

%%

subplot(3,1,1);
cycl=1:length(y{1});
hold on;
plot(cycl*3,y{1},"b","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf1,"r--","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf1_exo,"g--","LineWidth",2);
xlabel('Cycles');ylabel('Capacity');
legend('Exp','ARIMA','ARIMAX','location','best');
hold off;

subplot(3,1,2);
cycl=1:length(y{2});
hold on;
plot(cycl*3,y{2},"b","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf2,"r--","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf2_exo,"g--","LineWidth",2);
xlabel('Cycles');ylabel('Capacity');
% legend('Exp','ARIMA','ARIMAX');
hold off;

subplot(3,1,3)
cycl=1:length(y{3});
hold on;
plot(cycl*3,y{3},"b","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf3,"r--","LineWidth",2);
plot(cycl(Idx_est+1:end)*3,yf3_exo,"g--","LineWidth",2);
xlabel('Cycles');ylabel('Capacity');
% legend('Exp','ARIMA','ARIMAX');
hold off;

% figure();
% hold on
% plot(cycl,y1(1:length(cycl)),"b","LineWidth",2);
% plot(cycl,y2(1:length(cycl)),"b","LineWidth",2);
% plot(cycl,y3(1:length(cycl)),"b","LineWidth",2);
% plot(cycl(Idx_est+1:end),yf4,"r--","LineWidth",2);
