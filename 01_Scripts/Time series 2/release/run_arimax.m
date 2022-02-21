function [prediction]=run_arimax(runNo,pred_len,varargin)
%add path
    folder = fileparts(which(mfilename)); 
    addpath(genpath(folder));
    
    Br{1}=gdFun.Load_Multiple_Runs(runNo,false);
    n=1;% sampling freq of degradation data

    % get degradation trend cuve
    y(1)=Br{1}.cycleTable{2,'ahDchrge'};
        for j=2:height(Br{1}.cycleTable)/n      
            y(j)=Br{1}.cycleTable{n*j,'ahDchrge'};
             if y(j)<y(j-1)*0.7    % replace oddities in data where capacity suddenly drops for certain cycles
                y(j)= y(j-1);
             end
        end
    y(1)=Br{1}.cycleTable{2,'ahDchrge'};% first cycle capacity data is ususally wrong
    cell_cap=y(1);
    y=y./cell_cap;   
    yS=smooth(y,5,'lowess'); %smoothing
    
   if nargin<4
        train_len=length(y);
    elseif varargin{2}==0
        train_len=length(y);
    else
        train_len=floor(varargin{2}/n);
    end
         
    if nargin<3
        plot_test=false;
    else
        plot_test=varargin{1};
    end
    
    %convert to python data type and define model orders to train & compare
    X=py.numpy.array(yS(1:train_len));
    cycles=py.numpy.array(1:train_len);
    p_values = py.list([1:4]);
    d_values = py.list([1:2]);
    q_values = py.list([1:3]);
    pred_len=int16(pred_len/n);
    
    %find best model order based on rmse of prediction vs test data
    best_order=py.arima_functions.grid_search_orders(X,cycles,p_values,d_values,q_values);
    cycles=py.numpy.array(1:train_len+double(pred_len)); 
    %use the best model and all capacity data to make future predictions
    prediction=py.arima_functions.run_arima_prediction(X,cycles,pred_len,best_order);
    prediction=double(prediction);
    
%     if plot_pred==1
%         figure()
%         hold on;
%         plot(prediction(1,:)*n,prediction(2,:),'r--');
%         plot(1:n:n*length(y),y,'bs');
%         legend('forecast','data');
%         xlabel('Cycle');
%         ylabel('Normalised capacity');
%         hold off;
%     end
    
    if plot_test=='true'
        X2=py.numpy.array(yS);
        test=py.arima_functions.run_arima(X2,py.numpy.array(1:length(yS)),pred_len,best_order); 
        test=double(test);
        figure()
        hold on;
        plot(test(1,:)*n,test(2,:),'r--','LineWidth',1.5);
        plot(1:n:n*length(y),y,'bs');
        legend('Prediction on test data','data');
        xlabel('Cycle');
        ylabel('Normalised capacity');
        hold off;
    end
end