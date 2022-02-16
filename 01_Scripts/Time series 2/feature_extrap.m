function [Vsa_pred,Vmean_pred]=feature_extrap(train_len,pred_len,V_mean,V_sa)

     Vmean_train=py.numpy.array(V_mean(1:train_len));
     Vsa_train=py.numpy.array(V_sa(1:train_len));
     pred_len=int8(pred_len);
     best_order=[2,2,2]; % This should change with a grid search
     
    Vsa_pred=double(py.run_tseries.run_arima_prediction(Vsa_train,pred_len,best_order));
    Vmean_pred=double(py.run_tseries.run_arima_prediction(Vmean_train,pred_len,best_order));
    
end