%% ARIMAX model for capacity prediction

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

CapList={};
dVarQ={};
Cap_future={};
outlierCap={};
dNegPot={};
dChTime={};

%%
CapList={};
dVarQ={};
Cap_future={};
outlierCap={};
dChTime={};
n=3; %skip every n data points
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
    
    if i==4 | i==5 | i==6 | i==7 | i==8 | i==9
    y{i}(1:9)=[];
    end
    
   yS{i}=smooth(y{i},3); %smoothing
end

% 
% for i=1:9
% figure;
% hold on;
% % plot(yS{i},'b');
% plot(y{i},'r--');
% xlabel('Cycles');ylabel('%Capacity');
% hold off;
% end

%% ARIMA model

Idx_est=45;
Mdl=arima(2,1,0);


for i=1:length(Br)
[yP{i}]=polyfit(1:Idx_est,y{i}(1:Idx_est),1);
EstMdl(i)=estimate(Mdl,yS{i}(1:Idx_est));
end

for i=1:length(Br)
yf{i}=forecast(EstMdl(i),length(y{i})-Idx_est,yS{i}(1:Idx_est));
yP{i}=polyval(yP{i},1:length(y{i}));
cycl{i}=1:length(y{i});
end

subplot(4,2,1);
hold on;
plot(cycl{1}*n,y{1},"b","LineWidth",2);
plot(cycl{1}(Idx_est+1:end)*n,yf{1},"r--","LineWidth",2);
plot(cycl{1}*n,yP{1},"g--","LineWidth",2);
xlabel('Cycles');ylabel('%Capacity');
legend('Actual','ARIMA','Linear extrapolation','location','best');

for i=1:7
subplot(4,2,i+1);
hold on;
plot(cycl{3*i+1}*n,y{3*i+1},"b","LineWidth",2);
plot(cycl{3*i+1}(Idx_est+1:end)*n,yf{3*i+1},"r--","LineWidth",2);
plot(cycl{3*i+1}*n,yP{3*i+1},"g--","LineWidth",2);
xlabel('Cycles');ylabel('%Capacity');

end

hold off;
%%
%13 14 15 19 20 21 22 23 24
Mdl=arima(2,1,0);
plt=[13 14 15 19 20 21 22 23 24];

for i=1:length(plt)
 
Idx_est=length(y{plt(i)}); 
[yP{i}]=polyfit(1:Idx_est,y{plt(i)}(1:Idx_est),1);
EstMdl(i)=estimate(Mdl,yS{plt(i)}(1:Idx_est));    
yf{plt(i)}=forecast(EstMdl(i),50,y{plt(i)}(1:Idx_est));    
cycl{i}=1:length(y{plt(i)});  
cycl2{i}=cycl{i}(Idx_est):(cycl{i}(end)+49);    
yP{i}=polyval(yP{i},cycl2{i});
    
 
subplot(3,3, i)
hold on;
plot(cycl{i}*n,y{plt(i)},"b","LineWidth",2);
plot(cycl2{i}*n,yf{plt(i)},"r--","LineWidth",2);
plot(cycl2{i}*n,yP{i},"g--","LineWidth",2);

xlabel('Cycles');ylabel('%Capacity');

end
legend('Actual','ARIMA extrapolation');
hold off;

%%



% for i=1:3:length(Br)
% 
% hold on;
% plot(cycl{i}*n,y{i},"b","LineWidth",2);
% plot(cycl{i}(Idx_est+1:end)*n,yf{i},"r--","LineWidth",2);
% plot(cycl{i}*n,yP{i},"g--","LineWidth",2);
% xlabel('Cycles');ylabel('%Capacity');
% hold off;
% end

% figure;hold on;
% plot(cycl{2}*2,y{1},"b","LineWidth",2);
% plot(cycl{2}(Idx_est+1:end)*2,yf{2},"r--","LineWidth",2);
% xlabel('Cycles');ylabel('%Capacity');
% hold off;
