%% test deg exptrapolation based on ARIMA of features

% try with Elbrus data 
clear;
data2=gdFun.Load_Multiple_Runs([740,1018,1183,1190,1313,1554,1555],true);

% Extract discharge curves and 'normalised capacity' - stretch curves so that
% x axis goes from 0 to 1
% lines below extract all discharge curves and intepolate them to
% fixed length.If the data is already stitched then indexing will be
% simpler
%% Load capacity data
cell_cap=4.7;

for i=1:height(data2.RunData.cycleTable)    
    y(i)=data2.RunData.cycleTable{i,'ahDchrge'};
end

y=y./cell_cap;   
yS=y;
%% resample voltage curves 
Vtmp={};
Capcell={};
Vcell={};
Idx_disch=find((data2.RunData.dataTable{:,'currCell'}<0));
V_tmp=data2.RunData.dataTable{Idx_disch,'voltCell'};
Idx_bp=find(diff(V_tmp)>1);
%first disch curve
Vtmp{1}=V_tmp(1:Idx_bp(1)); Capcell{1}=(1:length(Vtmp{1}))./length(Vtmp{1});Vcell{1}=Vtmp{1};
%second to second-to-last
for i=2:length(Idx_bp)
    Vtmp{i}=V_tmp(Idx_bp(i-1)+1:Idx_bp(i));
    Captmp{i}=(1:length(Vtmp{i}))./length(Vtmp{i});
    Vcell{i}= Vtmp{i};
    Capcell{i}= Captmp{i};
end
%last discharge curve
Vtmp{i+1}=V_tmp(Idx_bp(i)+1:end);
Captmp{i+1}=(1:length(Vtmp{i+1}))./length(Vtmp{i+1}); 
Vcell{i+1}= Vtmp{i+1};
Capcell{i+1}= Captmp{i+1};   


 Vcell =Vcell(~cellfun('isempty',Vcell)) ; 
 Capcell =Capcell(~cellfun('isempty',Capcell)) ; 
 % resample so each section has the same length
 for i=1:length(Vcell)
     Vcell_intp{i}=interp1(Capcell{i},Vcell{i},linspace(Capcell{i}(1),1,500));
 end
Vcell=[Vcell_intp{:}];

clear ii Idx_disch V_tmp Vtmp Idx_bp Captmp Capcell Vcell_intp

%% seasonal decomposition and featurization

X=py.numpy.array(Vcell);
X_decomp=double(py.run_tseries.decomp_tseries(X));
% X_decomp(:,1:1000)=[]; %because trend value is alwasy wrong in first 1-2 cycle
Vseason=ones(500,length(X_decomp)/500);
for i=0:length(X_decomp)/500-1
    Vmean(i+1)=mean(X_decomp(1,500*i+1:500*(i+1)));    
    Vseason(:,i+1)=X_decomp(2,500*i+1:500*(i+1));
    Vsa(i+1)=trapz(abs(Vseason(:,i+1)-Vseason(:,1)))/500;
end

clear Vseason X X_decomp

%% Extrapolate features using arima
train_len=260;
Vmean_train=py.numpy.array(Vmean(1:train_len));
Vsa_train=py.numpy.array(Vsa(1:train_len));
pred_len=int8(150);
best_order=[2,2,2];
Vsa_pred=double(py.run_tseries.run_arima_prediction(Vsa_train,pred_len,best_order));
Vmean_pred=double(py.run_tseries.run_arima_prediction(Vmean_train,pred_len,best_order));

plot(Vsa);hold on;plot(Vsa_pred(1,:),Vsa_pred(2,:));hold off
figure
plot(Vmean);hold on;plot(Vmean_pred(1,:),Vmean_pred(2,:));hold off;


%% Extrapolate features using linear extrapolation

train_len=255;
Vmean_train=(Vmean(1:train_len));
Vsa_train=(Vsa(1:train_len));
fit_len=50;
pred_len=150;
x=1:train_len;
linear_fit_mean=polyfit(x(end-fit_len:end),Vmean_train(end-fit_len:end),1);
linear_fit_sa=polyfit(x(end-fit_len:end),Vsa_train(end-fit_len:end),1);


Vsa_pred=polyval(linear_fit_sa,length(x)+1:(length(x)+pred_len));
Vmean_pred=polyval(linear_fit_mean,length(x)+1:(length(x)+pred_len));

plot(Vsa);hold on;plot(length(x)+1:length(x)+pred_len,Vsa_pred(:));hold off
figure;
plot(Vmean);hold on;plot(length(x)+1:length(x)+pred_len,Vmean_pred(:)); hold off

%%  GPR
%normalise features
% [X_train,C,S]=normalize([double(Vmean_train)' double(Vsa_train)']);
X_train=[double(Vmean_train)' double(Vsa_train)'];
%GPR
gprMdl=fitrgp(X_train,yS(1:train_len)','Basis','linear','FitMethod','exact','PredictMethod','exact');

%Predict with training data
yhat = predict(gprMdl,X_train);
figure;
plot(yhat);hold on; plot(yS,'s');hold off;

%predict with extrapolation data and compare with real capacities
% X_val=normalize([Vmean_pred(2,:)' Vsa_pred(2,:)'],'center',C,'scale',S);


X_val=[Vmean_pred(2,:)' Vsa_pred(2,:)'];
yhat_v=predict(gprMdl,X_val);
figure
plot(train_len+1:train_len+double(pred_len),yhat_v,'--','LineWidth',2);hold on; plot(yS);
xlabel('Cycles');ylabel('Normalised capacity');legend('Prediction','Data');

%% Linear regression

X_train=[double(Vmean_train)' double(Vsa_train)'];
%lasso
[b1,FitInfo]=lasso(X_train,yS(1:train_len),'Alpha',0.7,'Intercept',true,'CV',10);
idxLambda1SE = FitInfo.Index1SE;
coef = b1(:,idxLambda1SE);
coef0 = FitInfo.Intercept(idxLambda1SE);
yhat = X_train*coef + coef0;

hold on;
figure;
plot(yhat,'--','LineWidth',2);hold on; plot(yS,'-');hold off;

% 
% X_val=[Vmean_pred(2,:)' Vsa_pred(2,:)'];
% yhat_v=X_val*coef + coef0;
% plot(train_len+1:train_len+double(pred_len),yhat_v);hold on; plot(yS);

X_val=[Vmean_pred' Vsa_pred'];
yhat_v=X_val*coef + coef0;
plot(train_len+1:train_len+pred_len,yhat_v);hold on; plot(yS);

%% random
% j=0;
% for i=1:length(Vcell)
%    j=j+length(Vcell{i});
% end
for i=1:height(data2.cycleTable)
    ismem(i)=ismember(i,data2.RunData.dataTable{:,'cycleNumber'});
    if ismem(i)==0
        i
    end
end


