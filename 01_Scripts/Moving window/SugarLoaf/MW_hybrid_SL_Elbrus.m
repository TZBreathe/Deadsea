%% Moving window GPR Elbrus data
%Load data
clear;
MaxCap=3;
runNumbers=[740 1018 1183 1190]; 
Br{1}=appendDegData(runNumbers);
runNumbers=[741 1019 1184 1191]; 
Br{2}=appendDegData(runNumbers);
runNumbers=[742 1020 1185 1192]; 
Br{3}=appendDegData(runNumbers);
runNumbers=[732 1013 1178 1187]; 
Br{4}=appendDegData(runNumbers);
runNumbers=[733 1014 1179 1188]; 
Br{5}=appendDegData(runNumbers);
runNumbers=[734 1015 1180 1189]; 
Br{6}=appendDegData(runNumbers);

runNumbers=[728 1012]; 
Br{7}=appendDegData(runNumbers);
runNumbers=[735 1016]; 
Br{8}=appendDegData(runNumbers);
runNumbers=[736 1017]; 
Br{9}=appendDegData(runNumbers);
runNumbers=[901 1033 1045 1196]; 
Br{10}=appendDegData(runNumbers);
runNumbers=[902 1034 1046 1197]; 
Br{11}=appendDegData(runNumbers);
runNumbers=[903 1035 1047 1198]; 
Br{12}=appendDegData(runNumbers);
runNumbers=[898 1030 1137]; 
Br{13}=appendDegData(runNumbers);
runNumbers=[899 1031 1138]; 
Br{14}=appendDegData(runNumbers);
runNumbers=[900 1032 1139]; 
Br{15}=appendDegData(runNumbers);
runNumbers=[895 1027 1042 1193]; 
Br{16}=appendDegData(runNumbers);
runNumbers=[896 1028 1043 1194]; 
Br{17}=appendDegData(runNumbers);

runNumbers=[897 1029 1044 1195]; 
Br{18}=appendDegData(runNumbers);
runNumbers=[1039 1143]; 
Br{19}=appendDegData(runNumbers);
runNumbers=[1040 1144]; 
Br{20}=appendDegData(runNumbers);
runNumbers=[1036 1140]; 
Br{21}=appendDegData(runNumbers);
runNumbers=[1037 1141]; 
Br{22}=appendDegData(runNumbers);
runNumbers=[1038 1142]; 
Br{23}=appendDegData(runNumbers);
runNumbers=[1041 1145]; 
Br{24}=appendDegData(runNumbers);

Br{25}=gdFun.Load_Breathe_Run(508);
Br{26}=gdFun.Load_Breathe_Run(509);




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
    CapList{i}(1)=CapList{i}(2);
end

Now=50;
Future=70;
Past=45;

%%
%var(dQ)
for i=1:length(Br)

    ii=1; %index for storing data points each cycle
    
    for j=Now:(length(CapList{i})-Future)  % index for 'now', begin at cycle 50 
    
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
 
% for i=1:length(Br)
%     ii=1;  
%     for j=Now:(length(CapList{i})-Now)    
%     k=1; % index 1-3 for storing avearge over 3 cycles
%         for l=(j-Past-2):(j-Past) %Average over 3 cycles for stable readings
%         idx_early=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==l]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
%             if isempty(idx_early) 
%              continue;
%             end
%          NegPot(k)=mean(Br{i}.RunData.dataTable{idx_early,'surpotNeg'});
%          k=k+1;
%         end
%      NegPotAvg_early=mean(NegPot);
%      
%     k=1; 
%         for l=(j-2):j  %Average over 3 cycles for stable readings
%          idx_now=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==l]&[Br{i}.RunData.dataTable{:,'currCell'}<0]);
%             if isempty(idx_now) 
%             continue;
%             end
%           NegPot(k)=mean(Br{i}.RunData.dataTable{idx_now,'surpotNeg'});
%             k=k+1;
%         end
%     NegPotAvg_now=mean(NegPot);
%     
%     dNegPot{i}(ii)=NegPotAvg_early-NegPotAvg_now;
%     ii=ii+1;
%     end
%     [dNegPot{i},outlierNegPot_tmp]=filloutliers(dNegPot{i},NaN); %set outliers to NaN to not affect regression
% %     outlierNegPot{i}=outlierNegPot_tmp;
% end
% 
% clear NegPotAvg_early NegPotAvg_now NegPot

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
idx=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
xdVarQ=[];xdNegPot=[];xdChTime=[];xCap_lag0=[];xCap_lag5=[];xCap_lag10=[];yCap=[];

%Arrange training data
for i =1:length(idx)
    xdVarQ=cat(1,xdVarQ,dVarQ{idx(i)}');
%     xdNegPot=cat(1,xdNegPot,dNegPot{idx(i)}');
    xdChTime=cat(1,xdChTime,dChTime{idx(i)}');
    xCap_lag0=cat(1,xCap_lag0,Cap_lag0{idx(i)}');
    xCap_lag5=cat(1,xCap_lag5,Cap_lag5{idx(i)}');
    xCap_lag10=cat(1,xCap_lag10,Cap_lag10{idx(i)}');
    yCap=cat(1,yCap,Cap_future{idx(i)}');
end

[X,C,S]=normalize([xdVarQ xdChTime xCap_lag0 xCap_lag5 xCap_lag10]);

%Fit GPR
gprMdl=fitrgp(X,yCap,'Basis','linear','FitMethod','exact','PredictMethod','exact');

%Predict with training data
yhat = predict(gprMdl,X);




%% Test regression

% define test data set
xdVarQ_v=[];xdChTime_v=[];xCap_lag0_v=[];xCap_lag5_v=[];xCap_lag10_v=[]; yCap_v=[];
idx_v=[25 26];
for i =1:length(idx_v)
    xdVarQ_v=cat(1,xdVarQ_v,dVarQ{idx_v(i)}');
%     xdNegPot=cat(1,xdNegPot,dNegPot{idx(i)}');
    xdChTime_v=cat(1,xdChTime_v,dChTime{idx_v(i)}');
    xCap_lag0_v=cat(1,xCap_lag0_v,Cap_lag0{idx_v(i)}'*MaxCap/4.8);
    xCap_lag5_v=cat(1,xCap_lag5_v,Cap_lag5{idx_v(i)}'*MaxCap/4.8);
    xCap_lag10_v=cat(1,xCap_lag10_v,Cap_lag10{idx_v(i)}'*MaxCap/4.8);
    yCap_v=cat(1, yCap_v,Cap_future{idx_v(i)}'*MaxCap/4.8);
end

X_v=normalize([xdVarQ_v xdChTime_v xCap_lag0_v xCap_lag5_v xCap_lag10_v],'center',C,'scale',S);

yhat_v = predict(gprMdl,X_v);


hold on;
scatter(yCap,yhat);
xlabel('True capacity');ylabel('Fitted capacity');
% hold off;
mean(abs(rmmissing(yhat-yCap)))
scatter(yCap_v,yhat_v);
xlabel('True capacity');ylabel('Fitted capacity');
plot(yCap_v,yCap_v,'red');
legend('Traning SL','Validation Elbrus','True Elbrus','True capacity','location','best');
hold off;
mean(abs(rmmissing(yhat_v-yCap_v)))




% 

