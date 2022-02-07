%% Run arima, linear, and power-law extrapolation of time-series degardation
% Set parameters
clear;
run_No=[1436];
pred_len=300;
fit_len=50;

%% make predictions
pred_arima=run_arima(run_No,pred_len,1);
[y,pred_l, pred_p]=run_extrap(run_No,fit_len,pred_len);

figure()
hold on;
plot(length(y):(length(y)+pred_len),pred_l(2,:),'r--');
% plot(length(y):(length(y)+pred_len),pred_p(2,:),'g--');
plot(length(y):3:(length(y)+pred_len-3),pred_arima(2,:),'bl--'); %mess
plot(1:length(y),y,'bs');
legend('linear extrapoltion','arima extrapolation','data');
xlabel('Cycle');
ylabel('Normalised capacity');
hold off;


%% temp for ibu
% clear;
% run_list=[1436 1437 1438 1439 1440 1441];
% pred_len=200;
% fit_len=50;
% figure();
% for i=1:length(run_list)
%     run_No=run_list(i);
% 
% 
% pred_arima=run_arima(run_No,pred_len);
% [y,pred_l, pred_p]=run_extrap(run_No,fit_len,pred_len);
% 
% subplot(2,3,i)
% hold on;
% plot(length(y):(length(y)+pred_len),pred_l(2,:),'r--');
% if i~=3 && i~=6
% plot(length(y):3:(length(y)+pred_len),pred_arima(2,:),'bl--'); %mess
% end
% plot(1:length(y),y,'bs');
% xlabel('Cycle');
% ylabel('Normalised capacity');
% if i==1
%     legend('linear extrapoltion','arima extrapolation','data','location','best');
% end
% title(run_No);
% ylim([0.95 1.02])
% hold off;
end

    