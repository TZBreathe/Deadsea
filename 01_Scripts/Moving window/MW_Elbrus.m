%% Extract features with 50 cycle difference and use as training set
% clear;
% 
% addpath(genpath('G:\09_Research_Development_Projects\00_Degradation data sets'));
% 
% runNumbers=[767 766 846]; 
% Br{1}=appendDegData(runNumbers);
% runNumbers=[765 764 845]; 
% Br{2}=appendDegData(runNumbers);
% runNumbers=[771 769 848]; 
% Br{3}=appendDegData(runNumbers);
% runNumbers=[776 773 850]; 
% Br{4}=appendDegData(runNumbers);
% 
% Br{5}=gdFun.Load_Breathe_Run(506);
% Br{6}=gdFun.Load_Breathe_Run(507);
% Br{7}=gdFun.Load_Breathe_Run(508);
% Br{8}=gdFun.Load_Breathe_Run(509);
% runNumbers=[808 807 852]; 
% Br{9}=appendDegData(runNumbers);
% runNumbers=[810 809 853]; 
% Br{10}=appendDegData(runNumbers);

MaxCap=4.9;

%%  Capacity list
% 
clear
MaxCap=4.9;
Br{1}=gdFun.Load_Breathe_Run(506);
Br{2}=gdFun.Load_Breathe_Run(507);
Br{3}=gdFun.Load_Breathe_Run(508);
Br{4}=gdFun.Load_Breathe_Run(509);
CapList={};
VarQ={};
Cap_future={};
dNegPot={};outlierNegPot={};
dChTime={};outlierdChTime={};
%%


for i=1:4
    
% create list of capacities to map to
CapList{i}(1)=Br{i}.RunData.cycleTable{1,'ahDchrge'};

    for j=2:height(Br{i}.RunData.cycleTable)      
        CapList{i}(j)=Br{i}.RunData.cycleTable{j,'ahDchrge'};
        if CapList{i}(j)<CapList{i}(j-1)/1.5    % replace oddities in data
            CapList{i}(j)= CapList{i}(j-1);
        end
    end


    %var(dQ)
    for j=50:1:(length(CapList{i})-50)            
    %Get indx for cycle j and cycle j-40 
    ind_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j-40]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
    ind_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&(Br{i}.RunData.dataTable{:,'currCell'}<0)); 
    
    %some datasets have missing sequences, fill these
    if length(ind_early)<10 | length(ind_now)<10
      VarQ{i}(j)=VarQ{i}(j-1);
      continue;
    end 
    
    VQ_early_tmp=[Br{i}.RunData.dataTable{ind_early,'voltCell'} Br{i}.RunData.dataTable{ind_early,'ahDchrge'}];
    [~,ia,ic]=unique(VQ_early_tmp(:,1));   %remove repeated V values
    VQ_early=VQ_early_tmp(ia,:);
    VQ_early(:,2)=VQ_early(:,2)./MaxCap; %Normalise capacity, make sense or not?
    %Interpolate to get same dimensions
    V_interp=linspace(VQ_early(1,1),VQ_early(end,1),1000);
    Q_interp=interp1(VQ_early(:,1),VQ_early(:,2),V_interp);
    VQ_early=[V_interp; Q_interp];
    
     VQ_now_tmp=[Br{i}.RunData.dataTable{ind_now,'voltCell'} Br{i}.RunData.dataTable{ind_now,'ahDchrge'}];
     [~,ia,ic]=unique(VQ_now_tmp(:,1));   %remove repeated V values
     VQ_now=VQ_now_tmp(ia,:);
     VQ_now(:,2)=VQ_now(:,2)./MaxCap;
     V_interp=linspace(VQ_now(1,1),VQ_now(end,1),1000);
     Q_interp=interp1(VQ_now(:,1),VQ_now(:,2),V_interp);
     VQ_now=[V_interp; Q_interp];
    
   %compute var(dQ)
     dQ=abs(VQ_early(2,:)-VQ_now(2,:));
     VarQ{i}(j)=var(dQ);
   Cap_future{i}(j)=CapList{i}(j+50)./MaxCap; % capacity at j+50
    end
    
    %Some cycle data is odd, replace these outliers, also get the index for
   %outliers

    [VarQ{i},outlier]=filloutliers(VarQ{i},median(VarQ{i}(50:end))); % replace outlier with median; should replace with zero when training
    Cap_future{i}(outlier)=median(Cap_future{i}(50:end)); %set a median capacity value for outliers, this should not be done if training model
end


% plot((VarQ(outlier)),Cap_future(outlier),'s');hold off;

   %% Avg negPot. Can use charge and/or discharge
for i=1:4
for j=50:1:(length(CapList{i})-50)    
    k=1;
    for l=(j-45):(j-40) %Average over 5 cycles for stable readings
      ind_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==l]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
       if isempty(ind_early) 
            continue;
      end
      NegPot(k)=mean(Br{i}.RunData.dataTable{ind_early,'surpotNeg'});
      k=k+1;
      
    end
     NegPotAvg_early=mean(NegPot);
     
    k=1; 
    for l=(j-5):j  %Average over 5 cycles for stable readings
      ind_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==l]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
       if isempty(ind_now) 
            continue;
      end
      NegPot(k)=mean(Br{i}.RunData.dataTable{ind_now,'surpotNeg'});
      k=k+1;
     end
    NegPotAvg_now=mean(NegPot);
    
    dNegPot{i}(j)=NegPotAvg_early-NegPotAvg_now;
    %Cap_future(j)=CapList(j+50)./MaxCap; % capacity at j+50
       
    
%     plot(NegPotAvg)
end

[dNegPot{i},outlierNegPottmp]=filloutliers(dNegPot{i},median(dNegPot{i}(50:end)));
% Cap_future{i}(outlierNegPot)=median(Cap_future{i}(50:end));
outlierNegPot{i}=outlierNegPottmp;


end


%%
% Totalcharge time
for i=1:4
dChTime{i}=ones(1,length(CapList)-50);  

for j=50:1:(length(CapList{i})-50)    

    
    ind_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j-40]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
    ind_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
    
      if isempty(ind_early) | isempty(ind_now)
       dChTime(j)=VarQ(j-1);
      continue;
      end 
    
    ChTime_early=Br{i}.RunData.dataTable{ind_early(end),'timeTest'}-Br{i}.RunData.dataTable{ind_early(1),'timeTest'};
    Chtime_now=Br{i}.RunData.dataTable{ind_now(end),'timeTest'}-Br{i}.RunData.dataTable{ind_now(1),'timeTest'};
    
    dChTime{i}(j)=Chtime_now/ChTime_early;
%     Cap_future(j)=CapList(j+50)./MaxCap; % capacity at j+50
    Cap_future{i}(j)=CapList{i}(j+50)./MaxCap; % capacity at j+50
    
end

[dChTime{i},outlierdChTime_tmp]=filloutliers(dChTime{i},median(dChTime{i}(50:end)));
outlierdChTime{i}=outlierdChTime_tmp;  


end
%% plot


figure()
plot(VarQ{1}(50:end),Cap_future{1}(50:end),'o');hold on;
plot(VarQ{2}(50:end),Cap_future{2}(50:end),'s');
plot(VarQ{3}(50:end),Cap_future{3}(50:end),'^');
plot(VarQ{4}(50:end),Cap_future{4}(50:end),'x');


figure();
plot(dNegPot{1}(~outlierNegPot{1}),Cap_future{1}(~outlierNegPot{1}),'o');hold on;
plot(dNegPot{2}(~outlierNegPot{2}),Cap_future{2}(~outlierNegPot{2}),'o');
plot(dNegPot{3}(~outlierNegPot{3}),Cap_future{3}(~outlierNegPot{3}),'o');
plot(dNegPot{4}(~outlierNegPot{4}),Cap_future{4}(~outlierNegPot{4}),'o');
hold off;
% 
figure();
plot(dChTime{1}(~outlierdChTime{1}),Cap_future{1}(~outlierdChTime{1}),'s');hold on;
plot(dChTime{2}(~outlierdChTime{2}),Cap_future{2}(~outlierdChTime{2}),'s');
plot(dChTime{3}(~outlierdChTime{3}),Cap_future{3}(~outlierdChTime{3}),'s');
plot(dChTime{4}(~outlierdChTime{4}),Cap_future{4}(~outlierdChTime{4}),'s');hold off;


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

