%% test deg exptrapolation based on ARIMA of features

% try with Elbrus data 
clear;
data=gdFun.Load_Breathe_Run( [898,1030,1137,1256,1528,1563]);
train_len=200;
pred_len=200;
n=1;

%% Load capacity data
i=0;
y=[];
cell_cap=2.8;
for j=1:length(data)
    
%     y(i+j)=data{j}.RunData.cycleTable{1,'ahDchrge'};
        for i=1:height(data{j}.RunData.cycleTable)  
            y=vertcat(y,data{j}.RunData.cycleTable{i,'ahDchrge'});
%             if y(i)>y(i-1)/1.05    % replace oddities in data where capacity suddenly drops for certain cycles
%                y(i)= [];
%             end
        end
end
     y=y./cell_cap;   
%      y(y(:)<0.5)=[];
%       y(51)=[];y(101-1)=[];y(151-2)=[];y(201-3)=[];y(251-4)=[];y(290-5)=[];y(340-6)=[];
%     %yS=smooth(y,5); %smoothing
      yS=y;
      
 %%
% Extract discharge curves and 'normalised capacity' - stretch curves so that
% x axis goes from 0 to 1
% lines below extract all discharge curves and intepolate them to
% fixed length.If the data is already stitched then indexing will be
% simpler
ii=0; %use for cycle index of Vcell
Vtmp={};
Vcell={};Capcell={};
for j=1:length(data)
   
    Idx_disch=find((data{1,j}.RunData.dataTable{:,'currCell'}<0));
    V_tmp=data{1,j}.RunData.dataTable{Idx_disch,'voltCell'};
    Idx_bp=find(diff(V_tmp)>0.8);
    %first disch curve
    Vtmp{1}=V_tmp(1:Idx_bp(1)); 
    ii=ii+1;
    Capcell{ii}=(1:length(Vtmp{1}))./length(Vtmp{1});Vcell{ii}=Vtmp{1};
    %second to second-to-last
    for i=2:length(Idx_bp)
        Vtmp{i}=V_tmp(Idx_bp(i-1)+1:Idx_bp(i));
        Captmp{i}=(1:length(Vtmp{i}))./length(Vtmp{i});
        ii=ii+1;
        Vcell{ii}= Vtmp{i};
        Capcell{ii}= Captmp{i};
%         if length(Vcell{ii})<1000 % remove partial disch
             %Vcell{ii}=[];
             %Capcell{ii}=[];
%             ii
%         end      
    end
    %last discharge curve
     Vtmp{i+1}=V_tmp(Idx_bp(i)+1:end);
     Captmp{i+1}=(1:length(Vtmp{i+1}))./length(Vtmp{i+1});
     ii=ii+1;
     Vcell{ii}= Vtmp{i+1};
     Capcell{ii}= Captmp{i+1};   
end

for i=1:length(Vcell)
    if length(Vcell{i})<1000
        Vcell{i}=[];
        Capcell{i}=[];
    end
end

 Vcell =Vcell(~cellfun('isempty',Vcell)) ; 
 Capcell =Capcell(~cellfun('isempty',Capcell)) ; 
 % resample so each section has the same length
 for i=1:length(Vcell)
     Vcell_intp{i}=interp1(Capcell{i},Vcell{i},linspace(Capcell{i}(1),1,500));
 end
Vcell=[Vcell_intp{:}];

% clear ii Idx_disch V_tmp Vtmp Idx_bp Captmp Capcell Vcell_intp

%% seasonal decomposition and featurization

X=py.numpy.array(Vcell);
X_decomp=double(py.arima_functions.decomp_tseries(X));
% X_decomp(:,1:1000)=[]; %because trend value is alwasy wrong in first 1-2 cycle
Vseason=ones(500,length(X_decomp)/500);
for i=0:length(X_decomp)/500-1
    Vmean(i+1)=mean(X_decomp(1,500*i+1:500*(i+1)));    
    Vseason(:,i+1)=X_decomp(2,500*i+1:500*(i+1));
    Vsa(i+1)=trapz(abs(Vseason(:,i+1)-Vseason(:,1)));
    Vsv(i+1)=var(abs(Vseason(:,i+1)-Vseason(:,1)))*100;
end

% clear Vseason X X_decomp

% save('temp35.mat')

%% Extrapolate features using arima
train_len=350;
pred_len=200;
n=3;
%yS=smooth(y,5); %smoothing
p_values = py.list([1:3]);
d_values = py.list([0:2]);
q_values = py.list([1:3]);

Vmean_train=log(Vmean(20:n:train_len));
Vsa_train=(Vsa(20:n:train_len));
Vsv_train=(Vsv(20:n:train_len));


Vmean_train=py.numpy.array(smooth(Vmean_train,8,'rloess'));
Vsa_train=py.numpy.array(smooth(Vsa_train,8,'rloess'));
Vsv_train=py.numpy.array(smooth(Vsv_train,8,'rloess')*100);
pred_len=int16(pred_len/n);

Vsa_pred=double(py.arima_functions.run_arima_auto(Vsa_train,pred_len));
Vmean_pred=double(py.arima_functions.run_arima_auto(Vmean_train,pred_len));
Vsv_pred=double(py.arima_functions.run_arima_auto(Vsv_train,pred_len));

plot(Vsa);hold on;plot(Vsa_pred(1,:),Vsa_pred(2,:));hold off
figure
plot(log(Vmean));hold on;plot(Vmean_pred(1,:),Vmean_pred(2,:));hold off;
figure
plot(Vsv*100);hold on;plot(Vsv_pred(1,:),Vsv_pred(2,:));hold off;

%% Extrapolate features using linear extrapolation
train_len=250;
pred_len=200;
pred_len=double(pred_len);
Vmean_train=(Vmean(1:train_len));
Vsa_train=(Vsa(1:train_len));
Vsv_train=(Vsv(1:train_len));


Vmean_train=smooth(Vmean_train,8,'rloess');
Vsa_train=smooth(Vsa_train,8,'rloess');
Vsv_train=smooth(Vsv_train,8,'rloess')*100;

fit_len=50;
x=1:train_len;
linear_fit_mean=polyfit(x(end-fit_len:end),Vmean_train(end-fit_len:end),1);
linear_fit_sa=polyfit(x(end-fit_len:end),Vsa_train(end-fit_len:end),1);
linear_fit_sv=polyfit(x(end-fit_len:end),Vsv_train(end-fit_len:end),1);


Vsa_pred=polyval(linear_fit_sa,length(x)+1:(length(x)+pred_len));
Vmean_pred=polyval(linear_fit_mean,length(x)+1:(length(x)+pred_len));
Vsv_pred=polyval(linear_fit_sv,length(x)+1:(length(x)+pred_len));

plot(Vsa);hold on;plot(length(x)+1:length(x)+pred_len,Vsa_pred(:));hold off
figure;
plot(Vmean);hold on;plot(length(x)+1:length(x)+pred_len,Vmean_pred(:)); hold off

%%  GPR
%normalise features
% [X_train,C,S]=normalize([double(Vmean_train)' double(Vsa_train)']);
X_train=[double(Vmean_train)' double(Vsa_train)' double(Vsv_train)'];
%GPR
gprMdl=fitrgp(X_train,yS(10:train_len)','Basis','linear','FitMethod','exact','PredictMethod','exact');

%Predict with training data
yhat = predict(gprMdl,X_train);
figure;
plot(yhat);hold on; plot(yS(1:1:end),'s');hold off;

%predict with extrapolation data and compare with real capacities
% X_val=normalize([Vmean_pred(2,:)' Vsa_pred(2,:)'],'center',C,'scale',S);
%  
% X_val=[Vmean_pred' Vsa_pred' Vsv_pred'];
% yhat_v=predict(gprMdl,X_val);
% figure
% plot(train_len+1:train_len+double(pred_len),yhat_v,'--','LineWidth',2);hold on; plot(yS);
% xlabel('Cycles');ylabel('Normalised capacity');legend('Prediction','Data');



X_val=[Vmean_pred(2,:)' Vsa_pred(2,:)' Vsv_pred(2,:)'];
yhat_v=predict(gprMdl,X_val);
figure
plot(train_len+1:train_len+double(pred_len),yhat_v,'--','LineWidth',2);hold on; plot(yS(1:1:end));
xlabel('Cycles');ylabel('Normalised capacity');legend('Prediction','Data');

%% Linear regression

yS(yS<0.6*yS(2))=[];
for i=2:length(yS)-1
    if yS(i+1)-yS(i)>3*(yS(i)-yS(i-1))
        yS(i+1)=yS(i);
    end
end

 X_train=[double(Vsa_train)' double(Vmean_train)'  double(Vsv_train)'];


%lasso
[b1,FitInfo]=lasso(X_train,yS(1:length(double(Vsa_train))),'Alpha',0.7,'Intercept',true,'CV',10);
idxLambda1SE = FitInfo.Index1SE;
coef = b1(:,idxLambda1SE);
coef0 = FitInfo.Intercept(idxLambda1SE);
yhat = X_train*coef + coef0;

hold on;

plot(yhat,'--','LineWidth',2);hold on; plot(yS,'-');hold off;

% X_val=[Vsa_pred(2,:)' Vmean_pred(2,:)'  Vsv_pred(2,:)'];
 X_val=[Vsa' log(Vmean)'  Vsv'*100];
yhat_v=X_val*coef + coef0;
plot(yhat_v);hold on; plot(yS);

% figure;
% X_val=[Vmean_pred' Vsa_pred' Vsv_pred'];
% yhat_v=X_val*coef + coef0;
% plot(train_len+1:train_len+pred_len,yhat_v);hold on; plot(yS);

