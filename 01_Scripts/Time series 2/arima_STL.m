function [yhat] = arima_STL(runNo,train_len,pred_len)
    % Get discharge curves and capacities
    [Vcell,yS]=get_disch_curves(runNo);
    
    % Get features from discharge curves
    [Vmean,Vsa]=get_features(Vcell);

    % Extrapolate features using arima
    [Vsa_pred,Vmean_pred]=feature_extrap(train_len,pred_len,Vmean,Vsa);

    % GPR
    Vmean_train=Vmean(1:train_len);
    Vsa_train=Vsa(1:train_len);
    X_train=[Vmean_train' Vsa_train'];
    gprMdl=fitrgp(X_train,yS(1:train_len)','Basis','linear','FitMethod','exact','PredictMethod','exact');

    %Prediction
    X_val=[Vmean_pred(2,:)' Vsa_pred(2,:)'];   
    yhat(:,2) = predict(gprMdl,X_val);
    yhat(:,1)= Vmean_pred(1,:)';
    
  
%     figure
%     plot(train_len+1:train_len+double(pred_len),yhat_v,'--','LineWidth',2);hold on; plot(yS);
%     xlabel('Cycles');ylabel('Normalised capacity');legend('Prediction','Data');




    % plot(Vsa);hold on;plot(Vsa_pred(1,:),Vsa_pred(2,:));hold off
    % figure
    % plot(Vmean);hold on;plot(Vmean_pred(1,:),Vmean_pred(2,:));hold off;

end



