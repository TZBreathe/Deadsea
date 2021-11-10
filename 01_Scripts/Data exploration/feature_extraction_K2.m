%% starter code for extracting features for degradation regression model

% K2 Br data
addpath(genpath('G:\09_Research_Development_Projects\00_Degradation data sets'));
filename=cell(6,1);
filename={'C46_FFC','C47_FFC','C48_lincc','C49_lincc','C51_lincc','C52_lincc'};
MaxCap=2.3;
MaxCurr=6.5;
%%

% variance of dQ(V)

for i=1:length(filename)
    
    Data=struct2cell(load(filename{i}));
    ind_10=find([Data{1, 1}.RunData.dataTable{:,'cycleNumber'}==7]&(Data{1, 1}.RunData.dataTable{:,'currCell'}<0));
    ind_50=find([Data{1, 1}.RunData.dataTable{:,'cycleNumber'}==52]&(Data{1, 1}.RunData.dataTable{:,'currCell'}<0)); %note differet cycles might be chosen here, 50,100,150, etc contain characterisation cycles and wont work
    
    VQ_10_tmp=[Data{1, 1}.RunData.dataTable{ind_10,'voltCell'} Data{1, 1}.RunData.dataTable{ind_10,'ahTotal'}]; 
    [~,ia,ic]=unique(VQ_10_tmp(:,1));   %remove repeated V values
    VQ_10=VQ_10_tmp(ia,:);
    VQ_10(:,2)=VQ_10(:,2)./MaxCap;
    
    V_interp=linspace(VQ_10(1,1),VQ_10(end,1),1000);
    Q_interp=interp1(VQ_10(:,1),VQ_10(:,2),V_interp);
    VQ_10=[V_interp; Q_interp];
    
    VQ_50_tmp=[Data{1, 1}.RunData.dataTable{ind_50,'voltCell'} Data{1, 1}.RunData.dataTable{ind_50,'ahTotal'}];
     [~,ia,ic]=unique(VQ_50_tmp(:,1));   %remove repeated V values
    VQ_50=VQ_50_tmp(ia,:);
    VQ_50(:,2)=VQ_50(:,2)./MaxCap;
     
    V_interp=linspace(VQ_50(1,1),VQ_50(end,1),1000);
    Q_interp=interp1(VQ_50(:,1),VQ_50(:,2),V_interp);
    VQ_50=[V_interp; Q_interp];
    
     dQ=abs(VQ_10(2,:)-VQ_50(2,:));
     VarQ(i)=var(dQ);
     
     EndCap(i)=Data{1, 1}.RunData.cycleTable{501,'ahDchrge'}./MaxCap; 
     
%         plot(VQ_10(1,:),VQ_10(2,:),VQ_50(1,:),VQ_50(2,:));
%        plot(dQ,VQ_10(1,:));

end
     

% plot(EndCap, VarQ,'O');

%%
% Average negPot


%  for i=3:length(filename)
%     
%     Data=struct2cell(load(filename{i}));
%     ind_10=find([Data{1, 1}.RunData.dataTable{:,'cycleNumber'}==7]&(Data{1, 1}.RunData.dataTable{:,'currCell'}>0));
%     ind_50=find([Data{1, 1}.RunData.dataTable{:,'cycleNumber'}==102]&(Data{1, 1}.RunData.dataTable{:,'currCell'}>0)); %note differet cycles might be chosen here, 50,100,150, etc contain characterisation cycles and wont work
%     
%     NegPotAvg_10=mean(Data{1, 1}.RunData.dataTable{ind_10,'surpotNeg'});
%     NegPotAvg_50=mean(Data{1, 1}.RunData.dataTable{ind_50,'surpotNeg'});
%     rNegPot(i)=NegPotAvg_50/NegPotAvg_10;
%     
%       EndCap(i)=Data{1, 1}.RunData.cycleTable{501,'ahDchrge'}./MaxCap;   
%  end
%  
%  
%  plot(rNegPot,EndCap);




% i=4; Data=struct2cell(load(filename{i}));
% j=1;k=2;idx(1)=1;
%  while j<height(Data{1, 1}.RunData.dataTable)
%      j=j+1;
%             
%      if Data{1, 1}.RunData.dataTable{j,'cycleNumber'}~=Data{1, 1}.RunData.dataTable{j-1,'cycleNumber'} 
%      idx(k)=j;
%      NegPotAvg(k)=mean(Data{1, 1}.RunData.dataTable{idx(k-1):idx(k),'surpotNeg'});
%      k=k+1;
%     
%      end
%       if k==50
%          break;
%      end
%  end
%  
 j=Data{1, 1}.RunData.dataTable{1,'cycleNumber'};
while j<50 
  ind=find([Data{1, 1}.RunData.dataTable{:,'cycleNumber'}==j]&(Data{1, 1}.RunData.dataTable{:,'currCell'}>0));
  NegPotAvg(j)=mean(Data{1, 1}.RunData.dataTable{ind,'surpotNeg'});
  j=j+1;
end


 NegPotAvg=filloutliers(NegPotAvg,'center'); 
 NegPotAvg=fillmissing(NegPotAvg,'next');

 NegPotAvg_10=mean(NegPotAvg(1:10));
 NegPotAvg_50=mean(NegPotAvg(end-10:end));
 dNegPot=NegPotAvg_10-NegPotAvg_50;

plot(NegPotAvg)

 
 %%
 % CC charge time
 
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

