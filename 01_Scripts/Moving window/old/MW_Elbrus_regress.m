%% Moving window prediction preliminary regression 

clear;
MaxCap=4.9;
% Br{1}=gdFun.Load_Breathe_Run(506);
% Br{2}=gdFun.Load_Breathe_Run(507);
% Br{3}=gdFun.Load_Breathe_Run(508);
% Br{4}=gdFun.Load_Breathe_Run(509);

runNumbers=[767 766 846]; 
Br{1}=appendDegData(runNumbers);
runNumbers=[765 764 845]; 
Br{2}=appendDegData(runNumbers);
runNumbers=[771 769 848]; 
Br{3}=appendDegData(runNumbers);
runNumbers=[776 773 850]; 
Br{4}=appendDegData(runNumbers);

Br{5}=gdFun.Load_Breathe_Run(506);
Br{6}=gdFun.Load_Breathe_Run(507);
Br{7}=gdFun.Load_Breathe_Run(508);
Br{8}=gdFun.Load_Breathe_Run(509);
CapList={};
VarQ={};outlierdQ={};
Cap_future={};outlierCap={};
dNegPot={};outlierNegPot={};
dChTime={};outlierdChTime={};


%% Extract features - List of capacities & Var(dQ)


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
%%
%var(dQ)
for i=1:length(Br)

    ii=1; %another index for storing data points each cycle
    
    for j=50:(length(CapList{i})-50)  % index for 'now', begin at cycle 50 
    
%     Cap_future{i}(ii)=CapList{i}(j+50); % Capacities of 'future' (j+50) cycles
      dCap{i}(ii)=CapList{i}(j)-CapList{i}(j+50); % Difference in capacity between now and 50 cycles later
    
        
    %Find early (j-40) and now (j) indexes for voltage and capacity
    ind_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j-40]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
    ind_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&(Br{i}.RunData.dataTable{:,'currCell'}<0)); 
    
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
    
    %Some cycle data contains oddities resulting in outliers in var(dQ),
    %find these outliers & their indexes
    %replace outlier with NaN so they don't affect regression. 
    [VarQ{i},outlierdQ_tmp]=filloutliers(VarQ{i},NaN); 
    outlierdQ{i}=outlierdQ_tmp;
     
end

   %% Avg negPot. Can use charge and/or discharge
   
 
for i=1:length(Br)
    ii=1;  
    for j=50:1:(length(CapList{i})-50)    
    k=1; % index 1-5 for storing avearge over 5 cycles
        for l=(j-44):(j-40) %Average over 5 cycles for stable readings
        ind_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==l]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
            if isempty(ind_early) 
             continue;
            end
         NegPot(k)=mean(Br{i}.RunData.dataTable{ind_early,'surpotNeg'});
         k=k+1;
        end
     NegPotAvg_early=mean(NegPot);
     
    k=1; 
        for l=(j-4):j  %Average over 5 cycles for stable readings
         ind_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==l]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
            if isempty(ind_now) 
            continue;
            end
          NegPot(k)=mean(Br{i}.RunData.dataTable{ind_now,'surpotNeg'});
            k=k+1;
        end
    NegPotAvg_now=mean(NegPot);
    
    dNegPot{i}(ii)=NegPotAvg_early-NegPotAvg_now;
    ii=ii+1;
    end
    [dNegPot{i},outlierNegPot_tmp]=filloutliers(dNegPot{i},NaN); %set outliers to NaN to not affect regression
    outlierNegPot{i}=outlierNegPot_tmp;
end


%%
% Totalcharge time
for i=1:length(Br)
 
    ii=1;
    for j=50:1:(length(CapList{i})-50)    

    ind_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j-40]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
    ind_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
    
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

    [dChTime{i},outlierdChTime_tmp]=filloutliers(dChTime{i},NaN);
    outlierdChTime{i}=outlierdChTime_tmp;  

end

%% check outliers

% Cap_outlier=find(outlierCap{1});
% dQ_outlier=find(outlierdQ{1});
% ChTime_outlier=find(outlierdChTime{1});
% NegPot_outlier=find(outlierNegPot{1});

%% multivariate regression with regularization
clear Br CapList ia ic 

xVarQ=normalize([VarQ{1}'; VarQ{2}'; VarQ{3}'; VarQ{4}';VarQ{5}';VarQ{6}';VarQ{7}';VarQ{8}']);
xdNegPot=normalize([dNegPot{1}';dNegPot{2}';dNegPot{3}';dNegPot{4}';dNegPot{5}';dNegPot{6}';dNegPot{7}';dNegPot{8}']);
xdChTime=normalize([dChTime{1}';dChTime{2}';dChTime{3}';dChTime{4}';dChTime{5}';dChTime{6}';dChTime{7}';dChTime{8}']);
% yCap_future=[Cap_future{1}'; Cap_future{2}'; Cap_future{3}'; Cap_future{4}';Cap_future{5}'; Cap_future{6}'; Cap_future{7}'; Cap_future{8}'];
ydCap=[dCap{1}'; dCap{2}'; dCap{3}'; dCap{4}';dCap{5}';dCap{6}'; dCap{7}'; dCap{8}'];
X=[ones(length(xVarQ),1) xVarQ xdNegPot xdChTime];

% lambda = 1e-04;
% b=lasso(X,yCap_future,'Lambda',lambda);
% b=mvregress(X,yCap_future);
b=mvregress(X,ydCap);
ycalc=X*b;

figure();
hold on;
% plot(yCap_future,ycalc,'o');hold on;
% plot(linspace(0,max(yCap_future),100),linspace(0,max(yCap_future),100));
plot(ydCap,ycalc,'o')
plot(linspace(0,max(ydCap),100),linspace(0,max(ydCap),100));
hold off;
% plot(X(:,3),yCap_future,'o');hold on;
% % plot(X(:,3),ycalc,'*'); hold off;
% plot(yCap_future,ycalc-yCap_future,'^');
mean(abs(rmmissing(ycalc-ydCap)))



%% Test regression
xVarQ=normalize([VarQ{5}'; VarQ{6}'; VarQ{7}'; VarQ{8}']);
xdNegPot=normalize([dNegPot{5}';dNegPot{6}';dNegPot{7}';dNegPot{8}']);
xdChTime=normalize([dChTime{5}';dChTime{6}';dChTime{7}';dChTime{8}']);
% yCap_future=[Cap_future{5}'; Cap_future{6}'; Cap_future{7}'; Cap_future{8}'];
ydCap=[dCap{5}';dCap{6}'; dCap{7}'; dCap{8}'];
X=[ones(length(xVarQ),1) xVarQ xdNegPot xdChTime];

% [b1,FitInfo]=lasso(X(:,2:4),yCap_future,'Alpha',0.7,'Intercept',true,'CV',10);
[b1,FitInfo]=lasso(X(:,2:4),ydCap,'Alpha',0.7,'Intercept',true,'CV',10);
idxLambda1SE = FitInfo.Index1SE;
coef = b1(:,idxLambda1SE);
coef0 = FitInfo.Intercept(idxLambda1SE);

yhat = X(:,2:4)*coef + coef0;
hold on;
% scatter(yCap_future,yhat);plot(yCap_future,yCap_future);
scatter(ydCap,yhat);plot(ydCap,yCap_future);
xlabel('True capacity');ylabel('Fitted capacity');
hold off;
mean(abs(rmmissing(yhat-yCap_future)))

% xVarQ=normalize([VarQ{1}';VarQ{2}';VarQ{3}';VarQ{4}']);
% yCap_future=([Cap_future{1}'; Cap_future{2}'; Cap_future{3}'; Cap_future{4}']);
% 
% X1=[ones(length(xVarQ),1) xVarQ];
% b1=mvregress(X1,yCap_future);
% ycalc1=X1*b1;
% 
% % b2=lasso(xVarQ,yCap_future,'Lambda',0.001,'Intercept',true);
% % ycalc2=X*b2;
% 
% % plot(xVarQ,yCap_future,'o');hold on;
% % plot(xVarQ,ycalc2,'.');
% % plot(xVarQ,ycalc,'--'); hold off;
% 
% plot(yCap_future,ycalc-yCap_future,'x');hold on;
% mean(abs(rmmissing(ycalc1-yCap_future)))




%% plot


% figure()
% plot((VarQ{1}),Cap_future{1},'o');hold on;
% plot((VarQ{2}),Cap_future{2},'s');
% plot(VarQ{3},Cap_future{3},'^');
% plot(VarQ{4},Cap_future{4},'x');
% % 
% % 
% figure();
% plot(dNegPot{1},Cap_future{1},'o');hold on;
% plot(dNegPot{2},Cap_future{2},'o');
% plot(dNegPot{3},Cap_future{3},'o');
% plot(dNegPot{4},Cap_future{4},'o');
% hold off;
% % 
% figure();
% plot(dChTime{1},Cap_future{1},'s');hold on;
% plot(dChTime{2},Cap_future{2},'s');
% plot(dChTime{3},Cap_future{3},'s');
% plot(dChTime{4},Cap_future{4},'s');hold off;


%% Avg disch T
% 
% for i=1:length(Br)
% 
%     k=1;
%     for j=10:15 %Average over cycle 10-15 for stable readings
%       idxT=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
%      
%      
%       Tavg(k)=mean(filloutliers(Br{i}.RunData.dataTable{idxT,'tempCell'},'center'));
%       k=k+1;
%       j=j+1;
%     end   
%    Tavg_10=mean(Tavg);
%    
%    k=1;
%    
%    for j=55:60 %Average over cycle 10-15 for stable readings
%       idxT=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
%       if isempty(ind_early) 
%       
%       continue;
%       end
%       Tavg(k)=mean(filloutliers(Br{i}.RunData.dataTable{idxT,'tempCell'},'center'));
%       k=k+1;
%       j=j+1;
%     end   
%     
%     Tavg_70=mean(Tavg);
%     
%     rTavg(i)=Tavg_10/Tavg_70;
%     dTavg(i)=Tavg_70-Tavg_10;
% end
% 
% plot(dTavg,EndCap,'s');
% 
% %See how T evolves
% % i=10;
% % j=1;
% % while j<149 
% %   ind=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
% %   Tavg(j)=mean(Br{i}.RunData.dataTable{ind,'tempCell'});
% %   j=j+1;
% % end
% % plot(Tavg);
% 
% 
% 

