%% starter code for extracting features for degradation regression model
% clear;
% Elbrus data more than 100 cycles
addpath(genpath('G:\09_Research_Development_Projects\00_Degradation data sets'));
%



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
runNumbers=[808 807 852]; 
Br{9}=appendDegData(runNumbers);
runNumbers=[810 809 853]; 
Br{10}=appendDegData(runNumbers);

MaxCap=4.9;

%% 
%variance of dQ(V)

for i=1:length(Br)
    
    
    ind_10=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==5]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
    ind_50=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==49]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
    
    VQ_10_tmp=[Br{i}.RunData.dataTable{ind_10,'voltCell'} Br{i}.RunData.dataTable{ind_10,'ahDchrge'}];
       [~,ia,ic]=unique(VQ_10_tmp(:,1));   %remove repeated V values
    VQ_10=VQ_10_tmp(ia,:);
    VQ_10(:,2)=VQ_10(:,2)./MaxCap; %Normalise capacity, make sense or not?
    
    V_interp=linspace(VQ_10(1,1),VQ_10(end,1),1000);
    Q_interp=interp1(VQ_10(:,1),VQ_10(:,2),V_interp);
    VQ_10=[V_interp; Q_interp];
    
    VQ_50_tmp=[Br{i}.RunData.dataTable{ind_50,'voltCell'} Br{i}.RunData.dataTable{ind_50,'ahDchrge'}];
     
     [~,ia,ic]=unique(VQ_50_tmp(:,1));   %remove repeated V values
     VQ_50=VQ_50_tmp(ia,:);
    VQ_50(:,2)=VQ_50(:,2)./MaxCap;
    
    V_interp=linspace(VQ_50(1,1),VQ_50(end,1),1000);
    Q_interp=interp1(VQ_50(:,1),VQ_50(:,2),V_interp);
    VQ_50=[V_interp; Q_interp];
    
     dQ=abs(VQ_10(2,:)-VQ_50(2,:));
     dQ= rmoutliers(dQ,'mean');
     VarQ(i)=rmoutliers(var(dQ),'mean');
     AvgQ(i)=rmoutliers(mean(dQ),'mean');
     MaxQ(i)=rmoutliers(max(dQ),'mean');
     
       EndCap(i)=Br{i}.RunData.cycleTable{149,'ahDchrge'}./MaxCap; 
     %     plot(VQ_10(1,:),VQ_10(2,:),VQ_50(1,:),VQ_50(2,:));
%        plot(dQ,VQ_50(1,:));

end
     

plot(EndCap, VarQ,'O');
% figure
% plot(EndCap, AvgQ,'O');
% figure
% plot(EndCap, MaxQ,'O');

%% Avg negPot. Can use charge and/or discharge

for i=1:length(Br)
    k=1;
    for j=10:15 %Average over cycle 10-15 for stable readings
      ind=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
      NegPot(k)=mean(Br{i}.RunData.dataTable{ind,'surpotNeg'});
      k=k+1;
      j=j+1;
    end
     NegPotAvg_10=mean(NegPot);
     k=1;
     
    for j=45:50  %Average over cycle 45-50 for stable readings
      ind=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&[Br{i}.RunData.dataTable{:,'currCell'}>0]);
      NegPot(k)=mean(Br{i}.RunData.dataTable{ind,'surpotNeg'});
      k=k+1;
      j=j+1;
    end
    NegPotAvg_50=mean(NegPot);
    dNegPot(i)=NegPotAvg_10-NegPotAvg_50;
     
%      NegPotAvg=filloutliers(NegPotAvg,'center'); 
%      NegPotAvg=fillmissing(NegPotAvg,'next');
%      NegPotAvg_10=mean(NegPotAvg(10:15));
%      NegPotAvg_50=mean(NegPotAvg(end-5:end));
%     plot(NegPotAvg)
end
plot(dNegPot,EndCap,'s');
%%
%CC charge time


for i=1:length(Br)

    MaxCurr_100=maxk(Br{i}.RunData.dataTable{:,'currCell'},100); %remove possible current spikes
    MaxCurr=MaxCurr_100(end);
    
    idx_10=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==10]&[Br{i}.RunData.dataTable{:,'currCell'}>MaxCurr*0.9]);
    idx_50=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==50]&[Br{i}.RunData.dataTable{:,'currCell'}>MaxCurr*0.9]);
    
    CCAh_10=Br{i}.RunData.dataTable{idx_10(end),'ahChrge'};
    CCAh_50=Br{i}.RunData.dataTable{idx_50(end),'ahChrge'};
    
    rCCAh(i)=CCAh_50/CCAh_10;
end

plot(rCCAh,EndCap,'s');

%% Avg disch T

for i=1:length(Br)

    k=1;
    for j=10:15 %Average over cycle 10-15 for stable readings
      idxT=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
      Tavg(k)=mean(filloutliers(Br{i}.RunData.dataTable{idxT,'tempCell'},'center'));
      k=k+1;
      j=j+1;
    end   
   Tavg_10=mean(Tavg);
   
   k=1;
   
   for j=55:60 %Average over cycle 10-15 for stable readings
      idxT=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
      Tavg(k)=mean(filloutliers(Br{i}.RunData.dataTable{idxT,'tempCell'},'center'));
      k=k+1;
      j=j+1;
    end   
    
    Tavg_70=mean(Tavg);
    
    rTavg(i)=Tavg_10/Tavg_70;
    dTavg(i)=Tavg_70-Tavg_10;
end

plot(dTavg,EndCap,'s');

%See how T evolves
% i=10;
% j=1;
% while j<149 
%   ind=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==j]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
%   Tavg(j)=mean(Br{i}.RunData.dataTable{ind,'tempCell'});
%   j=j+1;
% end
% plot(Tavg);




