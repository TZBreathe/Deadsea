clear
i=1;
pred_len=120;
fit_len=50;
train_len=100;
run_No=[740,1018,1183,1190,1313,1554,1555];

Br=gdFun.Load_Multiple_Runs(run_No,false);


%%
pred_len=120;
n=1;
pred_len=int16(pred_len/n);
y(1)=Br.cycleTable{2,'ahDchrge'};
        for j=2:height(Br.cycleTable)/n      
            y(j)=Br.cycleTable{n*j,'ahDchrge'};
             if y(j)<y(j-1)/1.5    % replace oddities in data where capacity suddenly drops for certain cycles
              y(j)= y(j-1);
             end
        end
y(1)=Br.cycleTable{2,'ahDchrge'};% first cycle capacity data is ususally wrong
yS=smooth(y,1); %smoothing
    %convert to python data type and define model orders to train & compare
cycles=15:length(yS)';
cycle=py.numpy.array(cycles);
X=py.numpy.array(yS(15:end));

p_values = py.list([1:4]);
d_values = py.list([0:2]);
q_values = py.list([1:3]);
best_order=py.arima_functions.grid_search_orders(X,cycle,p_values,d_values,q_values);


pred_arima=py.arima_functions.run_arima_ex(X,cycle,pred_len,best_order);
pred_arima=double(pred_arima);


figure()
hold on;
plot(pred_arima(1,:),pred_arima(2,:),'k--'); %mess
plot(1:length(yS),yS,'bs');
legend('arima','data');
%%
pred_len=120;
n=3;
pred_len=int16(pred_len/n);
y(1)=Br.cycleTable{2,'ahDchrge'};
        for j=2:height(Br.cycleTable)/n      
            y(j)=Br.cycleTable{n*j,'ahDchrge'};
             if y(j)<y(j-1)/1.5    % replace oddities in data where capacity suddenly drops for certain cycles
              y(j)= y(j-1);
             end
        end
y(1)=Br.cycleTable{2,'ahDchrge'};% first cycle capacity data is ususally wrong
yS=smooth(y,5); %smoothing
    %convert to python data type and define model orders to train & compare
cycles=1:length(yS)';
cycle=py.numpy.array(cycles);
X=py.numpy.array(yS);
p_values = py.list([1:4]);
d_values = py.list([0:2]);
q_values = py.list([1:3]);
best_order=py.arima_functions.grid_search_orders(X,cycle,p_values,d_values,q_values);

cycles=[1:length(yS)+double(pred_len)]';
cycle=py.numpy.array(cycles);
pred_arima=py.arima_functions.run_arima_prediction(X,cycle,pred_len,best_order);
pred_arima=double(pred_arima);
figure()
 hold on;
plot(pred_arima(1,:),pred_arima(2,:),'k--'); %mess
plot(1:length(yS),yS,'bs');
legend('arima','arima_{ex}','data');