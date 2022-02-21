function [yhat,Vmean,Vsa,Vsv,IR] = run_STL(runNo,pred_len,varargin)
%add path
    folder = fileparts(which(mfilename)); 
    addpath(genpath(folder));    

% Get discharge curves and capacities (and IR drop at beginning of
% discharge)

    [Vcell,yS,IR]=get_disch_curves(runNo);
    
% check inputs    
    
    if nargin<4
        train_len=length(yS);
    elseif varargin{2}==0 || varargin{2}>length(yS)
        train_len=length(yS);
    else
        train_len=floor(varargin{2}/n);
    end
         
    if nargin<3
        plot_test=false;
    else
        plot_test=varargin{1};
    end
    
% Get 3 features from discharge curves
    [Vmean,Vsa,Vsv]=get_features(Vcell);
    
% Extrapolate features
    %sometimes feature sample size != capacity sample size due to
    %inconsistencies in data stitching or outiler handlling. To avoid bug
    %the smaller of the two are used in training models
    train_len=min(length(Vmean),length(yS)); 
    [Vmean_pred,Vsa_pred,Vsv_pred,IR_pred]=feature_extrap(Vmean,Vsa,Vsv,IR,pred_len,train_len);

    % GPR training
    Vmean_train=Vmean(1:train_len);
    Vsa_train=Vsa(1:train_len);
    Vsv_train=Vsv(1:train_len);
    IR_train=IR(1:train_len)';
    X_train=[Vmean_train Vsa_train Vsv_train IR_train];
    gprMdl=fitrgp(X_train,yS(1:train_len),'Basis','linear','FitMethod','exact','PredictMethod','exact');
    %ytrain=predict(gprMdl,X_train);
    
    % GPR Prediction
    X_val=[Vmean_pred' Vsa_pred' Vsv_pred' IR_pred'];   
    yhat(:,2) = predict(gprMdl,X_val);
    yhat(:,1)= [train_len+1:train_len+pred_len]';
    
     try
         if plot_test=='true'
            train_len2=floor(train_len*0.75);
            pred_len2=ceil(train_len*0.25);
            X_train2=X_train(1:train_len2,:);
            [Vmean_pred2,Vsa_pred2,Vsv_pred2,IR_pred2]=feature_extrap(Vmean,Vsa,Vsv,IR,pred_len2,train_len2);
            X_val2=[Vmean_pred2' Vsa_pred2' Vsv_pred2' IR_pred2'];  
            gprMdl2=fitrgp(X_train2,yS(1: train_len2),'Basis','linear','FitMethod','exact','PredictMethod','exact');
            test(:,2)=predict(gprMdl,X_val2);
            test(:,1)= [train_len2+1:train_len2+pred_len2]';

            figure()
            hold on;
            plot(test(:,1),test(:,2),'r--','LineWidth',1.5);
            plot(1:n:n*length(yS),yS,'bs');
            legend('Prediction on test data','data');
            xlabel('Cycle');
            ylabel('Normalised capacity');
            hold off;
         end
     catch
     end
    
end



