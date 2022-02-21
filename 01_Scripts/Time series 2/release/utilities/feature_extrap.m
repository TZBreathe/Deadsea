function [Vmean_pred,Vsa_pred,Vsv_pred,IR_pred]=feature_extrap(Vmean,Vsa,Vsv,IR,pred_len,train_len)

%ARIMA prediction - not robust 
%     pred_len=int16(pred_len);
%     Vmean_train=py.numpy.array(smooth(Vmean(1:train_len),8,'rlowess'));
%     Vsa_train=py.numpy.array(smooth(Vsa(1:train_len),8,'rlowess'));
%     Vsv_train=py.numpy.array(smooth(Vsv(1:train_len),8,'rlowess')*10);
%     IR_train=py.numpy.array(smooth(IR(1:train_len),8,'rlowess')*100);
% 
%     Vmean_pred=double(py.arima_functions.run_arima_auto(Vmean_train,pred_len));
%     Vsa_pred=double(py.arima_functions.run_arima_auto(Vsa_train,pred_len));
%     Vsv_pred=double(py.arima_functions.run_arima_auto(Vsv_train,pred_len));
%     IR_pred=double(py.arima_functions.run_arima_auto(IR_train,pred_len));

%linear extrapolation - more robust but only linear trend possible
      fit_len=75;
      Vmean_train=Vmean(1:train_len);
      Vsa_train=Vsa(1:train_len);
      Vsv_train=Vsv(1:train_len);
      IR=smooth(IR,'rlowess');
      IR_train=IR(1:train_len);
      
      x=1:train_len;
      linear_fit_mean=polyfit(x(end-fit_len:end),Vmean_train(end-fit_len:end),1);
      linear_fit_sa=polyfit(x(end-fit_len:end),Vsa_train(end-fit_len:end),1);
      linear_fit_sv=polyfit(x(end-fit_len:end),Vsv_train(end-fit_len:end),1);
      linear_fit_ir=polyfit(x(end-fit_len:end),IR_train(end-fit_len:end),1);

      Vsa_pred=polyval(linear_fit_sa,length(x)+1:(length(x)+pred_len));
      Vmean_pred=polyval(linear_fit_mean,length(x)+1:(length(x)+pred_len));
      Vsv_pred=polyval(linear_fit_sv,length(x)+1:(length(x)+pred_len));  
      IR_pred=polyval(linear_fit_ir,length(x)+1:(length(x)+pred_len));  
    
%     figure  
%     plot(Vsa);hold on;plot(length(x)+1:length(x)+pred_len,Vsa_pred(:));hold off
%     figure;
%     plot(Vmean);hold on;plot(length(x)+1:length(x)+pred_len,Vmean_pred(:)); hold off
%     figure
%     plot(IR);hold on;plot(length(x)+1:length(x)+pred_len,IR_pred(:)); hold off
      
end